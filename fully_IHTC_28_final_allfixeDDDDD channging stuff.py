
import gurobipy as gp
from gurobipy import GRB
import json
from collections import defaultdict

# === Load data ===
with open("i03.json") as f:
    data = json.load(f)

# === Setup ===
days = data["days"]
days_list = list(range(days))
shift_types = data["shift_types"]
shifts_list = list(range(len(shift_types)))
age_order = {ag: i for i, ag in enumerate(data["age_groups"])}
weights = data["weights"]

# === Extract core entities ===
patients = data["patients"]
occupants = data["occupants"]
rooms = data["rooms"]
ots = data["operating_theaters"]
surgeons = data["surgeons"]
nurses = data["nurses"]

patient_ids = [p["id"] for p in patients]
room_ids = [r["id"] for r in rooms]
ot_ids = [ot["id"] for ot in ots]
surgeon_ids = [s["id"] for s in surgeons]
nurse_ids = [n["id"] for n in nurses]

# === Parse patient, room, nurse, surgeon data ===
room_capacity = {r["id"]: r["capacity"] for r in rooms}
ot_avail = {ot["id"]: ot["availability"] for ot in ots}
surgeon_avail = {s["id"]: s["max_surgery_time"] for s in surgeons}
nurse_skill = {n["id"]: n["skill_level"] for n in nurses}
nurse_working = defaultdict(set)
nurse_maxload = {}
for n in nurses:
    for ws in n["working_shifts"]:
        nurse_working[n["id"]].add((ws["day"], shift_types.index(ws["shift"])))
        nurse_maxload[(n["id"], ws["day"], shift_types.index(ws["shift"]))] = ws["max_load"]

# Patient data
p_release = {p["id"]: p["surgery_release_day"] for p in patients}
p_due = {p["id"]: p.get("surgery_due_day", days - 1) for p in patients}
p_los = {p["id"]: p["length_of_stay"] for p in patients}
p_mand = {p["id"]: p.get("mandatory", False) for p in patients}
p_gender = {p["id"]: p["gender"] for p in patients}
p_age = {p["id"]: p["age_group"] for p in patients}
p_surgdur = {p["id"]: p["surgery_duration"] for p in patients}
p_surgeon = {p["id"]: p["surgeon_id"] for p in patients}
p_work = {p["id"]: p["workload_produced"] for p in patients}
p_skill = {p["id"]: p["skill_level_required"] for p in patients}
p_incomp_rooms = {p["id"]: set(p.get("incompatible_room_ids", [])) for p in patients}

# Occupant data
occ_gender = {o["id"]: o["gender"] for o in occupants}
occ_age = {o["id"]: o["age_group"] for o in occupants}
occ_work = {o["id"]: o["workload_produced"] for o in occupants}
occ_skill = {o["id"]: o["skill_level_required"] for o in occupants}
occ_room_day = defaultdict(list)
for o in occupants:
    for d in range(o["length_of_stay"]):
        occ_room_day[(o["room_id"], d)].append(o["id"])

# === Build Model ===
model = gp.Model("IHTC2024_Revised")
model.setParam("OutputFlag", 1)

# === Variables ===
admit = model.addVars(patient_ids, days_list, vtype=GRB.BINARY, name="admit")
scheduled = model.addVars(patient_ids, vtype=GRB.BINARY, name="scheduled")
room_assign = model.addVars(patient_ids, room_ids, vtype=GRB.BINARY, name="room_assign")
theater_assign = model.addVars(patient_ids, ot_ids, days_list, vtype=GRB.BINARY, name="theater_assign")
nurse_assign = model.addVars(nurse_ids, room_ids, days_list, shifts_list, vtype=GRB.BINARY, name="nurse_assign")

# === Admission constraints ===
for p in patient_ids:
    release, due = p_release[p], p_due[p]
    model.addConstr(gp.quicksum(admit[p, d] for d in range(release, due + 1)) == scheduled[p])
    if p_mand[p]:
        model.addConstr(scheduled[p] == 1)

# === Room assignment ===
for p in patient_ids:
    model.addConstr(gp.quicksum(room_assign[p, r] for r in room_ids) == scheduled[p])
    for r in p_incomp_rooms[p]:
        model.addConstr(room_assign[p, r] == 0)

# === Prevent nurse assignment to unavailable shifts ===
for n in nurse_ids:
    for d in days_list:
        for s in shifts_list:
            if (d, s) not in nurse_working[n]:
                for r in room_ids:
                    model.addConstr(nurse_assign[n, r, d, s] == 0)


# === Helper function to create linearized AND variable ===
def linearize_and(model, x, y, name):
    z = model.addVar(vtype=GRB.BINARY, name=name)
    model.addConstr(z <= x)
    model.addConstr(z <= y)
    model.addConstr(z >= x + y - 1)
    return z

# === H3: Surgeon daily time limit ===
for s in surgeon_ids:
    for d in days_list:
        model.addConstr(
            gp.quicksum(
                p_surgdur[p] * theater_assign[p, t, d]
                for p in patient_ids if p_surgeon[p] == s
                for t in ot_ids
            ) <= surgeon_avail[s][d]
        )

# === H4: OT daily time availability ===
for t in ot_ids:
    for d in days_list:
        model.addConstr(
            gp.quicksum(p_surgdur[p] * theater_assign[p, t, d] for p in patient_ids)
            <= ot_avail[t][d]
        )

# === H7: Room capacity constraint ===
for r in room_ids:
    for d in days_list:
        model.addConstr(
            gp.quicksum(
                room_assign[p, r] * admit[p, a]
                for p in patient_ids
                for a in range(p_release[p], p_due[p] + 1)
                if a <= d < a + p_los[p]
            ) + len(occ_room_day.get((r, d), [])) <= room_capacity[r]
        )

# === A4: Surgery must happen on admission day ===
for p in patient_ids:
    for d in days_list:
        model.addConstr(gp.quicksum(theater_assign[p, t, d] for t in ot_ids) == admit[p, d])

# === A5: Room must be covered by at least one nurse if there's any workload ===
for r in room_ids:
    for d in days_list:
        for s in shifts_list:
            workload_terms = []
            workload_exists = False
            for p in patient_ids:
                los = p_los[p]
                for a in range(p_release[p], p_due[p] + 1):
                    if a <= d < a + los:
                        shift_idx = 3 * (d - a) + s
                        workload_terms.append(p_work[p][shift_idx] * admit[p, a] * room_assign[p, r])
                        workload_exists = True
            for occ in occ_room_day.get((r, d), []):
                shift_idx = 3 * d + s
                workload_terms.append(occ_work[occ][shift_idx])
                workload_exists = True
            if workload_exists:
                model.addConstr(
                    gp.quicksum(nurse_assign[n, r, d, s] for n in nurse_ids) >= 1,
                    name=f"nurse_cover_{r}_{d}_{s}"
                )

# === A3: Prevent assigning nurses to rooms with zero workload ===
for n in nurse_ids:
    for r in room_ids:
        for d in days_list:
            for s in shifts_list:
                workload_terms = []
                for p in patient_ids:
                    los = p_los[p]
                    for a in range(p_release[p], p_due[p] + 1):
                        if a <= d < a + los:
                            shift_idx = 3 * (d - a) + s
                            workload_terms.append(p_work[p][shift_idx] * admit[p, a] * room_assign[p, r])
                for occ in occ_room_day.get((r, d), []):
                    shift_idx = 3 * d + s
                    workload_terms.append(occ_work[occ][shift_idx])
                if workload_terms:
                    workload_indicator = model.addVar(vtype=GRB.BINARY)
                    model.addConstr(workload_indicator >= gp.quicksum(workload_terms) / 1000)
                    model.addConstr(nurse_assign[n, r, d, s] <= workload_indicator)

# === H1: No Gender Mixing (Hard Constraint) ===
gender_in_room = model.addVars(room_ids, days_list, ["A", "B"], vtype=GRB.BINARY, name="gender_in_room")

for r in room_ids:
    for d in days_list:
        # Patients assigned
        for p in patient_ids:
            los = p_los[p]
            gender = p_gender[p]
            for a in range(p_release[p], p_due[p] + 1):
                if a <= d < a + los:
                    and_var = linearize_and(model, admit[p, a], room_assign[p, r], f"gender_{p}_{r}_{d}")
                    model.addConstr(gender_in_room[r, d, gender] >= and_var)

        # Occupants
        for occ in occ_room_day.get((r, d), []):
            gender = occ_gender[occ]
            gender_in_room[r, d, gender].LB = 1  # occupant imposes presence of gender

        # Enforce no mix
        model.addConstr(gender_in_room[r, d, "A"] + gender_in_room[r, d, "B"] <= 1)
# === Soft Constraints & Objective ===
age_penalties = []
skill_mismatch_penalties = []
overload_penalties = []
ot_open_vars = {}
surgeon_transfer_penalties = []
admission_delays = []



# === S1: Age Group Mixing Penalty ===
for r in room_ids:
    for d in days_list:
        age_vars = []
        for p in patient_ids:
            los = p_los[p]
            for a in range(p_release[p], p_due[p] + 1):
                if a <= d < a + los:
                    and_var = linearize_and(model, admit[p, a], room_assign[p, r], f"agmix_{p}_{r}_{d}")
                    age_vars.append((age_order[p_age[p]], and_var))
        for occ in occ_room_day.get((r, d), []):
            age_vars.append((age_order[occ_age[occ]], 1))
        if age_vars:
            max_age = model.addVar(lb=0, ub=2, vtype=GRB.INTEGER)
            min_age = model.addVar(lb=0, ub=2, vtype=GRB.INTEGER)
            model.addGenConstrMax(max_age, [ag for ag, _ in age_vars])
            model.addGenConstrMin(min_age, [ag for ag, _ in age_vars])
            age_penalties.append(max_age - min_age)



# === S2: Skill mismatch penalty (correct linearization) ===
for n in nurse_ids:
    skill = nurse_skill[n]
    for d in days_list:
        for s in shifts_list:
            for r in room_ids:
                for p in patient_ids:
                    los = p_los[p]
                    for a in range(p_release[p], p_due[p] + 1):
                        if a <= d < a + los:
                            index = 3 * (d - a) + s
                            req_skill = p_skill[p][index]
                            diff = max(0, req_skill - skill)
                            if diff > 0:
                                and1 = linearize_and(model, admit[p, a], room_assign[p, r], f"mismatch1_{p}_{r}_{d}_{s}")
                                and2 = linearize_and(model, and1, nurse_assign[n, r, d, s], f"mismatch2_{p}_{r}_{n}_{d}_{s}")
                                skill_mismatch_penalties.append(diff * and2)
                for occ in occ_room_day.get((r, d), []):
                    index = 3 * d + s
                    req_skill = occ_skill[occ][index]
                    diff = max(0, req_skill - skill)
                    if diff > 0:
                        skill_mismatch_penalties.append(diff * nurse_assign[n, r, d, s])
for n in nurse_ids:
    skill = nurse_skill[n]
    for d in days_list:
        for s in shifts_list:
            for r in room_ids:
                for p in patient_ids:
                    los = p_los[p]
                    for a in range(p_release[p], p_due[p] + 1):
                        if a <= d < a + los:
                            index = 3 * (d - a) + s
                            req_skill = p_skill[p][index]
                            diff = max(0, req_skill - skill)
                            if diff > 0:
                                and1 = linearize_and(model, admit[p, a], room_assign[p, r], f"mismatch1_{p}_{r}_{d}_{s}")
                                and2 = linearize_and(model, and1, nurse_assign[n, r, d, s], f"mismatch2_{p}_{r}_{n}_{d}_{s}")
                                skill_mismatch_penalties.append(diff * and2)
                for occ in occ_room_day.get((r, d), []):
                    index = 3 * d + s
                    req_skill = occ_skill[occ][index]
                    diff = max(0, req_skill - skill)
                    if diff > 0:
                        skill_mismatch_penalties.append(diff * nurse_assign[n, r, d, s])
for n in nurse_ids:
    skill = nurse_skill[n]
    for d in days_list:
        for s in shifts_list:
            for r in room_ids:
                for p in patient_ids:
                    los = p_los[p]
                    for a in range(p_release[p], p_due[p] + 1):
                        if a <= d < a + los:
                            index = 3 * (d - a) + s
                            req_skill = p_skill[p][index]
                            diff = max(0, req_skill - skill)
                            if diff > 0:
                                and1 = linearize_and(model, admit[p, a], room_assign[p, r], f"mismatch1_{p}_{r}_{d}_{s}")
                                and2 = linearize_and(model, and1, nurse_assign[n, r, d, s], f"mismatch2_{p}_{r}_{n}_{d}_{s}")
                                skill_mismatch_penalties.append(diff * and2)
                for occ in occ_room_day.get((r, d), []):
                    index = 3 * d + s
                    req_skill = occ_skill[occ][index]
                    diff = max(0, req_skill - skill)
                    if diff > 0:
                        pass
# === S4: Nurse overload  ===
for n in nurse_ids:
    for d in days_list:
        for s in shifts_list:
            workload_terms = []
            for r in room_ids:
                room_workload = 0
                for p in patient_ids:
                    los = p_los[p]
                    for a in range(p_release[p], p_due[p] + 1):
                        if a <= d < a + los:
                            shift_idx = 3 * (d - a) + s
                            room_workload += p_work[p][shift_idx]
                for occ in occ_room_day.get((r, d), []):
                    shift_idx = 3 * d + s
                    room_workload += occ_work[occ][shift_idx]
                workload_terms.append(room_workload * nurse_assign[n, r, d, s])
            if (n, d, s) in nurse_maxload:
                maxload = nurse_maxload[(n, d, s)]
                total_workload = model.addVar(lb=0.0, name=f"total_workload_{n}_{d}_{s}")
                model.addConstr(total_workload == gp.quicksum(workload_terms))
                over = model.addVar(lb=0.0, name=f"overload_{n}_{d}_{s}")
                model.addConstr(over >= total_workload - maxload)
                overload_penalties.append(over)

# === S5: Open OT penalty ===
for t in ot_ids:
    for d in days_list:
        open_var = model.addVar(vtype=GRB.BINARY, name=f"openOT_{t}_{d}")
        surgeries = gp.quicksum(theater_assign[p, t, d] for p in patient_ids)
        model.addConstr(open_var >= surgeries / 1000)
        ot_open_vars[t, d] = open_var

# === S6: Surgeon transfer penalty ===
for s in surgeon_ids:
    for d in days_list:
        ot_used = []
        for t in ot_ids:
            used = model.addVar(vtype=GRB.BINARY, name=f"surgeonUsed_{s}_{t}_{d}")
            assigned = [theater_assign[p, t, d] for p in patient_ids if p_surgeon[p] == s]
            if assigned:
                model.addConstr(used >= gp.quicksum(assigned) / 1000)
            ot_used.append(used)
        surgeon_transfer_penalties.append(gp.quicksum(ot_used) - 1)

# === S7: Admission delay ===
for p in patient_ids:
    if not p_mand[p]:
        delay = gp.quicksum(d * admit[p, d] for d in range(p_release[p], days))
        admission_delays.append(delay - p_release[p] * scheduled[p])


# === S8: Optional patients are not scheduled ===
unscheduled = [1 - scheduled[p] for p in patient_ids if not p_mand[p]]

# === S2: Continuity of Care === not using S2, as looping through everything is not efficient
nurse_cares = model.addVars(patient_ids, nurse_ids, vtype=GRB.BINARY, name="nurse_cares")
continuity_penalties = []
for p in patient_ids:
    for n in nurse_ids:
        for r in room_ids:
            for d in days_list:
                for s in shifts_list:
                    for a in range(p_release[p], p_due[p] + 1):
                        if a <= d < a + p_los[p]:
                            and1 = linearize_and(model, admit[p, a], room_assign[p, r], f"cont_and1_{p}_{r}_{a}")
                            and2 = linearize_and(model, and1, nurse_assign[n, r, d, s], f"cont_and2_{p}_{n}_{r}_{d}_{s}")
                            model.addConstr(nurse_cares[p, n] >= and2)
    # Penalty: total number of distinct nurses assigned to patient p
    continuity_penalties.append(gp.quicksum(nurse_cares[p, n] for n in nurse_ids))


# === Objective  ===
model.setObjective(
    weights["room_mixed_age"] * gp.quicksum(age_penalties) +
    weights["room_nurse_skill"] * gp.quicksum(skill_mismatch_penalties) +
    weights["nurse_eccessive_workload"] * gp.quicksum(overload_penalties) +
    weights["open_operating_theater"] * gp.quicksum(ot_open_vars.values()) +
    weights["surgeon_transfer"] * gp.quicksum(surgeon_transfer_penalties) +
    weights["patient_delay"] * gp.quicksum(admission_delays) +
    weights["unscheduled_optional"] * gp.quicksum(unscheduled) +
    weights["continuity_of_care"] * gp.quicksum(continuity_penalties)
    ,
    GRB.MINIMIZE
)

model.setParam("TimeLimit", 600)
model.setParam("Presolve", 2)
model.setParam("LogFile", "gurobi_log.txt")
model.optimize() #optimizing the model here

if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
    final_patients = []
    for p in patient_ids:
        for d in days_list:
            if admit[p, d].X > 0.5:
                assigned_room = next((r for r in room_ids if room_assign[p, r].X > 0.5), None)
                assigned_ot = next((t for t in ot_ids if theater_assign[p, t, d].X > 0.5), None)
                if assigned_room and assigned_ot:
                    final_patients.append({
                        "id": p,
                        "admission_day": d,
                        "room": assigned_room,
                        "operating_theater": assigned_ot
                    })

    final_nurses = []
    for n in nurse_ids:
        nurse_assignments = []
        for d in days_list:
            for s in shifts_list:
                # Always include the shift, even if no rooms are assigned
                assigned_rooms = [r for r in room_ids if nurse_assign[n, r, d, s].X > 0.5]
                nurse_assignments.append({
                    "day": d,
                    "shift": shift_types[s],
                    "rooms": assigned_rooms  # include empty list ok
                })
        final_nurses.append({
            "id": n,
            "assignments": nurse_assignments
        })

    final_solution = {
        "patients": final_patients,
        "nurses": final_nurses
    }

    with open("sol_i12.json", "w") as f:
        json.dump(final_solution, f, indent=2)

    print("✅ Final solution exported to 'sol_i12.json'")
else:
    print("❌ No feasible solution found.")
#---- printing the summary from here ----
print("\n📝 Summary of Assignments")

print("\n👨‍⚕️ Patients:")
for entry in final_patients:
    print(f"  - Patient {entry['id']} admitted on day {entry['admission_day']} to room {entry['room']} in OT {entry['operating_theater']}")

print("\n👩‍⚕️ Nurse Assignments:")
for nurse in final_nurses:
    print(f"  - Nurse {nurse['id']}:")
    for a in nurse["assignments"]:
        day = a["day"]
        shift = a["shift"]
        rooms = ", ".join(a["rooms"]) if a["rooms"] else "(no rooms)"
        print(f"    - Day {day}, Shift {shift}: {rooms}")
