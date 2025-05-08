enter here a description on what is what and what is happening where

important: need gurobi academic license

# IHTC 2024 – Integrated Healthcare Timetabling

This repository contains my submission for the Integrated Healthcare Timetabling Competition (IHTC 2024). It includes a Mixed Integer Programming (MIP) model and Constraint Programming (CP) models for solving instances of the integrated healthcare scheduling problem as defined in the [official problem description](https://ihtc2024.github.io).

## Contents

- `fully_IHTC_28_final_allfixeDDDDD.py`: MIP model implemented using Gurobi.
- `IHTC_CP_FINAL_no_heuristics.py`: CP model using Google OR-Tools (basic).
- `IHTC_CP_FINAL_heuristics.py`: CP model with added search heuristics.
- `i03.json`: Example instance file.
- `sol_i12.json` / `solution.json`: Output examples.
- `ihtc2024_problem_specification.pdf`: Official problem spec.

## Requirements

### MIP (Gurobi)
- Python 3.8+
- `gurobipy`
- Valid Gurobi academic license

### CP (OR-Tools)
- Python 3.8+
- `ortools`:
  ```bash
  pip install ortools
