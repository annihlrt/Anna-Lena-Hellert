
# Seminar paper IHTC 2024

This repository contains my code. It includes a Mixed Integer Programming MIP model and Constraint Programming CP models for solving instances of the integrated healthcare scheduling problem. Here more info and all the data: (https://ihtc2024.github.io).

## Contents

- `IHTC_MIP_final.py`: MIP model implemented using Gurobi.
- `IHTC_CP_FINAL_no_heuristics.py`: CP model using Google OR-Tools (basic).
- `IHTC_CP_FINAL_heuristics.py`: CP model with added search heuristics.
- `i03.json`: Example instance file.
- `IHTP_Validator.exe`and `validator.cc` for validation. can also be found on IHTC2024 website

## Requirements

### MIP (Gurobi)
- Python 3.8+
- `gurobipy`
- Valid Gurobi academic license

### CP (OR-Tools)
- Python 3.8+
- `ortools`

### Validation
- eg. GNU compiler g++

## Instances
- For my paper only the included instances are necessary. But for completion all instances were tested.
- Additional public instances (i01–i30) can be downloaded from the [official competition page](https://ihtc2024.github.io).

### Validation
To validate results the following code needs to be written into the terminal. 
````bash
./IHTP_Validator instance.json solution.json verbose
