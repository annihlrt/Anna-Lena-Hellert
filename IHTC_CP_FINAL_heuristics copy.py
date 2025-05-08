
from ortools.sat.python import cp_model
import json

# Load the JSON input
with open("i02.json") as f:
    data = json.load(f)

model = cp_model.CpModel()

D = data["days"]
S = len(data["shift_types"])
day_range = range(D)
shift_range = range(S)
patients = data["patients"]
nurses = data["nurses"]
surgeons = data["surgeons"]
ots = data["operating_theaters"]
rooms = data["rooms"]
occupants = data["occupants"]
room_ids = [r["id"] for r in rooms]
ot_ids = [t["id"] for t in ots]
nurse_ids = [str(n["id"]) for n in nurses]
patient_ids = [p["id"] for p in patients]
shifts = data["shift_types"]
room_cap = {r['id']: r['capacity'] for r in rooms}

# Normalize nurse IDs to strings
for nurse in nurses:
    nurse["id"] = str(nurse["id"])

# Occupants in room/day map
occ_map = {}
for occ in occupants:
    for d in range(occ["length_of_stay"]):
        occ_map.setdefault((occ["room_id"], d), []).append(occ)

# Variables
admit = {(p["id"], d): model.NewBoolVar(f"admit_{p['id']}_{d}") for p in patients for d in day_range}
scheduled = {p["id"]: model.NewBoolVar(f"scheduled_{p['id']}") for p in patients}
room_assign = {(p["id"], r): model.NewBoolVar(f"room_{p['id']}_{r}") for p in patients for r in room_ids}
surgery_assign = {(p["id"], t, d): model.NewBoolVar(f"surg_{p['id']}_{t}_{d}") for p in patients for t in ot_ids for d in day_range}

# Nurse assign only for defined shifts
working_shifts = set((n["id"], ws["day"], shifts.index(ws["shift"])) for n in nurses for ws in n["working_shifts"])
nurse_assign = {
    (n["id"], r["id"], d, s): model.NewBoolVar(f"nurse_{n['id']}_{r['id']}_{d}_{s}")
    for n in nurses
    for r in rooms
    for d in day_range
    for s in shift_range
    if (n["id"], d, s) in working_shifts
}

# Patient presence tracking
# Correct patient presence tracking
presence = {}
for p in patients:
    pid = p["id"]
    LOS = p["length_of_stay"]
    release = p["surgery_release_day"]
    due = p.get("surgery_due_day", D - 1)

    for r in room_ids:
        for d in day_range:
            var = model.NewBoolVar(f"present_{pid}_{r}_{d}")
            presence[(pid, r, d)] = var

            # Compute admissible admission day exactly LOS days before or earlier
            possible_admit_days = [a for a in range(release, min(due, d) + 1) if a <= d < a + LOS]

            # Patient is present in room r on day d exactly if they were admitted on exactly one of these days and assigned to room r.
            admits_and_room = [model.NewBoolVar(f"admits_and_room_{pid}_{r}_{d}_{a}") for a in possible_admit_days]
            
            for idx, a in enumerate(possible_admit_days):
                model.AddBoolAnd([admit[pid, a], room_assign[pid, r]]).OnlyEnforceIf(admits_and_room[idx])
                model.AddBoolOr([admit[pid, a].Not(), room_assign[pid, r].Not()]).OnlyEnforceIf(admits_and_room[idx].Not())
            
            # Patient presence equals the OR over these admission-day-room combinations
            if admits_and_room:
                model.AddBoolOr(admits_and_room).OnlyEnforceIf(var)
                for admit_room_var in admits_and_room:
                    model.AddImplication(admit_room_var, var)
            else:
                model.Add(var == 0)



# === Constraints ===


# H5 + H6: Mandatory patients must be scheduled within release/due window
for p in patients:
    pid = p["id"]
    release = p["surgery_release_day"]
    due = p.get("surgery_due_day", D - 1)
    model.Add(sum(admit[pid, d] for d in range(release, due + 1)) == scheduled[pid])
    if p.get("mandatory", False):
        model.Add(scheduled[pid] == 1)
        # Hx: Exactly one room assigned if patient is scheduled
for p in patients:
    pid = p["id"]
    model.Add(sum(room_assign[pid, r] for r in room_ids) == scheduled[pid])


# H2: Incompatible rooms
for p in patients:
    for r in p.get("incompatible_room_ids", []):
        model.Add(room_assign[p["id"], r] == 0)

# H3: Surgeon daily limits
for s in surgeons:
    sid = s["id"]
    for d in day_range:
        model.Add(sum(
            p["surgery_duration"] * surgery_assign[p["id"], t, d]
            for p in patients if p["surgeon_id"] == sid for t in ot_ids
        ) <= s["max_surgery_time"][d])

# H4: OT daily capacity
for t in ot_ids:
    ot = next(ot for ot in ots if ot["id"] == t)
    for d in day_range:
        model.Add(sum(
            p["surgery_duration"] * surgery_assign[p["id"], t, d] for p in patients
        ) <= ot["availability"][d])

# A4: Surgery must be scheduled on admission day
for p in patients:
    pid = p["id"]
    for d in day_range:
        model.Add(sum(surgery_assign[pid, t, d] for t in ot_ids) == admit[pid, d])

# === H7: Room capacity (REVISED) ===
for r in room_ids:
    for d in day_range:
        occ = len(occ_map.get((r, d), []))
        patient_in_room = [presence[pid, r, d] for pid in patient_ids if (pid, r, d) in presence]
        model.Add(sum(patient_in_room) + occ <= room_cap[r])

# === H1: No gender mix (SIMPLIFIED) ===
for r in room_ids:
    for d in day_range:
        A_present = model.NewBoolVar(f"has_A_{r}_{d}")
        B_present = model.NewBoolVar(f"has_B_{r}_{d}")

        A_sources = []
        B_sources = []

        for p in patients:
            pid = p["id"]
            if (pid, r, d) in presence:
                if p["gender"] == "A":
                    A_sources.append(presence[pid, r, d])
                elif p["gender"] == "B":
                    B_sources.append(presence[pid, r, d])

        for occ in occ_map.get((r, d), []):
            if occ["gender"] == "A":
                A_sources.append(model.NewConstant(1))
            elif occ["gender"] == "B":
                B_sources.append(model.NewConstant(1))

        if A_sources:
            model.AddMaxEquality(A_present, A_sources)
        else:
            model.Add(A_present == 0)

        if B_sources:
            model.AddMaxEquality(B_present, B_sources)
        else:
            model.Add(B_present == 0)

        # Cannot have both A and B present
        model.Add(A_present + B_present <= 1)

# === A5: Assign nurses when patients are present (REVISED) ===
for r in room_ids:
    for d in day_range:
        for s in shift_range:
            needs_nurse = []
            for p in patients:
                pid = p["id"]
                if (pid, r, d) in presence:
                    needs_nurse.append(presence[pid, r, d])
            for occ in occ_map.get((r, d), []):
                needs_nurse.append(model.NewConstant(1))
            if needs_nurse:
                required = model.NewBoolVar(f"nurse_required_{r}_{d}_{s}")
                model.AddMaxEquality(required, needs_nurse)
                available = [nurse_assign[n, r, d, s] for n in nurse_ids if (n, r, d, s) in nurse_assign]
                if available:
                    model.Add(sum(available) >= 1).OnlyEnforceIf(required)

# === Objective ===
objective_terms = []
weights = data["weights"]
age_order = {ag: i for i, ag in enumerate(data["age_groups"])}
# === Soft Constraints Continued ===
# S1: Age mix in rooms
for r in room_ids:
    for d in day_range:
        ages = []
        for p in patients:
            pid = p["id"]
            if (pid, r, d) in presence:
                ages.append(model.NewIntVar(0, len(age_order), f"age_{pid}_{r}_{d}"))
        for occ in occ_map.get((r, d), []):
            ages.append(age_order[occ["age_group"]])
        if ages:
            max_age = model.NewIntVar(0, len(age_order), f"max_age_{r}_{d}")
            min_age = model.NewIntVar(0, len(age_order), f"min_age_{r}_{d}")
            model.AddMaxEquality(max_age, ages)
            model.AddMinEquality(min_age, ages)
            diff = model.NewIntVar(0, len(age_order), f"age_diff_{r}_{d}")
            model.Add(diff == max_age - min_age)
            objective_terms.append(diff * weights["room_mixed_age"])


# S2: Nurse skill mismatch
nurse_skill = {n["id"]: n["skill_level"] for n in nurses}
for n in nurse_ids:
    skill = nurse_skill[n]
    for r in room_ids:
        for d in day_range:
            for s in shift_range:
                if (n, r, d, s) in nurse_assign:
                    penalty_terms = []
                    for p in patients:
                        pid = p["id"]
                        for a in range(p["surgery_release_day"], p.get("surgery_due_day", D - 1) + 1):
                            if a <= d < a + p["length_of_stay"]:
                                shift_idx = 3 * (d - a) + s
                                if shift_idx < len(p["skill_level_required"]):
                                    req = p["skill_level_required"][shift_idx]
                                    if req > skill:
                                        active = model.NewBoolVar(f"active_skill_{pid}_{r}_{d}_{s}_{n}")
                                        model.AddBoolAnd([admit[pid, a], room_assign[pid, r], nurse_assign[n, r, d, s]]).OnlyEnforceIf(active)
                                        penalty_terms.append(active * (req - skill))
                    if penalty_terms:
                        mismatch = model.NewIntVar(0, 1000, f"mismatch_{n}_{r}_{d}_{s}")
                        model.Add(mismatch == sum(penalty_terms))
                        objective_terms.append(mismatch * weights["room_nurse_skill"])

# S3: Continuity of care
for p in patients:
    pid = p["id"]
    nurse_present = []
    for n in nurse_ids:
        contributes = []
        for r in room_ids:
            for d in day_range:
                for s in shift_range:
                    for a in range(p["surgery_release_day"], p.get("surgery_due_day", D - 1) + 1):
                        if a <= d < a + p["length_of_stay"]:
                            if (n, r, d, s) in nurse_assign:
                                shift_idx = 3 * (d - a) + s
                                if shift_idx < len(p["workload_produced"]) and p["workload_produced"][shift_idx] > 0:
                                    contributes.append(nurse_assign[n, r, d, s])
        if contributes:
            contrib = model.NewBoolVar(f"nurse_{n}_used_for_{pid}")
            model.AddMaxEquality(contrib, contributes)
            nurse_present.append(contrib)
    if nurse_present:
        total_nurses = model.NewIntVar(0, len(nurse_present), f"nurse_count_{pid}")
        model.Add(total_nurses == sum(nurse_present))
        objective_terms.append(total_nurses * weights["continuity_of_care"])



# S8: Optional patients not admitted
for p in patients:
    if not p.get("mandatory", False):
        pid = p["id"]
        not_admitted = model.NewBoolVar(f"not_admitted_{pid}")
        model.Add(not_admitted == 1 - scheduled[pid])
        objective_terms.append(not_admitted * weights["unscheduled_optional"])

# S7: Admission delay
#for p in patients:
 #   pid = p["id"]
 #   release = p["surgery_release_day"]
 #   delay = model.NewIntVar(0, D, f"delay_{pid}")
 #   model.Add(delay == sum((d - release) * admit[pid, d] for d in range(release, D)))
 #   objective_terms.append(delay * weights["patient_delay"])

for p in patients:
    pid = p["id"]
    release = p["surgery_release_day"]
    for d in range(release, D):
        if (pid, d) in admit:
            delay = model.NewIntVar(0, D, f"delay_{pid}")
            model.Add(delay == d - release).OnlyEnforceIf(admit[pid, d])
            model.Add(delay == 0).OnlyEnforceIf(admit[pid, d].Not())
            objective_terms.append(delay * weights["patient_delay"])
# S4: Nurse workload excess
# Default nurse workload limit per shift
default_nurse_shift_limit = 10  # adjust based on data scale

nurse_limits = {}
for n in nurse_ids:
    nurse_limits[n] = {}
    for d in range(D):  # days
        for s in range(S):  # shift types
            nurse_limits[n][(d, s)] = default_nurse_shift_limit
# Construct limits: (nurse_id, day, shift) → max load
for n in nurse_ids:
    for d in day_range:
        for s in shift_range:
            if (n, d, s) in nurse_limits[n]:
                load_terms = []
                for r in room_ids:
                    if (n, r, d, s) in nurse_assign:
                        for p in patients:
                            pid = p["id"]
                            for a in range(p["surgery_release_day"], min(D, p.get("surgery_due_day", D - 1) + 1)):
                                if a <= d < a + p["length_of_stay"]:
                                    shift_idx = 3 * (d - a) + s
                                    if shift_idx < len(p["workload_produced"]):
                                        load = p["workload_produced"][shift_idx]
                                        if load > 0:
                                            active = model.NewBoolVar(f"active_load_{pid}_{n}_{r}_{d}_{s}")
                                            model.AddBoolAnd([
                                                admit[pid, a],
                                                room_assign[pid, r],
                                                nurse_assign[n, r, d, s]
                                            ]).OnlyEnforceIf(active)
                                            load_terms.append(active * load)

                if load_terms:
                    total_load = model.NewIntVar(0, 1000, f"total_load_{n}_{d}_{s}")
                    model.Add(total_load == sum(load_terms))
                    max_allowed = nurse_limits[n][(d, s)]
                    excess = model.NewIntVar(0, 1000, f"excess_{n}_{d}_{s}")
                    model.Add(excess == total_load - max_allowed).OnlyEnforceIf(total_load > max_allowed)
                    model.Add(excess == 0).OnlyEnforceIf(total_load <= max_allowed)
                    capped = model.NewIntVar(0, 5, f"capped_excess_{n}_{d}_{s}")
                    model.AddMinEquality(capped, [excess, model.NewConstant(5)])
                    objective_terms.append(capped * weights["nurse_eccessive_workload"])

# === Heursitics ===
# Heuristic: Prioritize tight-window (inflexible) patients to be admitted earlier


# Step 1: Create (pid, d, flexibility) tuples
admit_vars = []
for p in patients:
    pid = p["id"]
    release = p["surgery_release_day"]
    due = p.get("surgery_due_day", D - 1)
    flexibility = due - release
    for d in range(release, due + 1):
        admit_vars.append((flexibility, pid, d, admit[pid, d]))

# Step 2: Sort by (flexibility ASC, day ASC)
admit_vars.sort(key=lambda x: (x[0], x[2]))

# Step 3: Extract variables only
ordered_admit_vars = [var for _, _, _, var in admit_vars]

# Step 4: Apply decision strategy
model.AddDecisionStrategy(
    ordered_admit_vars,
    cp_model.CHOOSE_FIRST,
    cp_model.SELECT_MIN_VALUE
)

# Heuristic: Prioritize assigning high-skilled nurses first
nurse_skill = {n["id"]: n["skill_level"] for n in nurses}

# Build a list of (skill, nurse_id, room_id, day, shift, BoolVar)
skill_ordered_assignments = []
for n in nurse_ids:
    skill = nurse_skill[n]  # assume this is a precomputed dict
    for r in sorted(room_ids):
        for d in range(D):
            for s in range(S):
                key = (n, r, d, s)
                if key in nurse_assign:
                    skill_ordered_assignments.append(
                        (skill, n, r, d, s, nurse_assign[key])
                    )

# Sort by skill DESC, then by nurse, room, day, shift
skill_ordered_assignments.sort(key=lambda x: (-x[0], x[1], x[2], x[3], x[4]))

# Extract only the BoolVars
ordered_nurse_vars = [var for *_, var in skill_ordered_assignments]

# Apply decision strategy
model.AddDecisionStrategy(
    ordered_nurse_vars,
    cp_model.CHOOSE_FIRST,
    cp_model.SELECT_MAX_VALUE
)

# === Final Objective ===
model.Minimize(sum(objective_terms))

# === Solver and JSON Output ===
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = 300
solver.parameters.log_search_progress = True  # Print intermediate search logs  
solver.parameters.use_lns = True
solver.parameters.num_search_workers = 4  # or more depending on your machine
status = solver.Solve(model)

if status == cp_model.OPTIMAL:
    print("Optimal solution found.")
elif status == cp_model.FEASIBLE:
    print("Feasible solution found (not necessarily optimal).")
else:
    print("No solution found.")

if status in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
    patient_output = []
    for p in patients:
        pid = p["id"]
        if solver.Value(scheduled[pid]):
            admission_day = next(d for d in day_range if solver.Value(admit[pid, d]))
            room = next(r for r in room_ids if solver.Value(room_assign[pid, r]))
            ot = next(t for t in ot_ids for d in day_range if solver.Value(surgery_assign[pid, t, d]))
            patient_output.append({
                "id": pid,
                "admission_day": admission_day,
                "room": room,
                "operating_theater": ot
            })

    nurse_output = []
    for n in nurse_ids:
        assignments = []
        for d in day_range:
            for s in shift_range:
                rooms = [r for r in room_ids if (n, r, d, s) in nurse_assign and solver.Value(nurse_assign[n, r, d, s])]
                assignments.append({
                    "day": d,
                    "shift": shifts[s],
                    "rooms": rooms
                })
        nurse_output.append({"id": n, "assignments": assignments})

    output = {"patients": patient_output, "nurses": nurse_output}

    with open("solution.json", "w") as f:
        json.dump(output, f, indent=2)
    print("Solution written to solution.json")
else:
    print("No solution found.")