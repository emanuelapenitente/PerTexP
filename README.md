# PerTexP (Pertussis Time Exploration) — MATLAB codes
PerTexP is a MATLAB tool allowing the scenario-based exploration of pertussis dynamics under maternal and infant vaccinations
# PerTexP 

This repository contains MATLAB scripts and functions to **simulate different vaccination scenarios** and to **compute the control reproduction number** \(\mathcal{R}_c\) under different vaccination strategies.

The implementation follows the model described in the accompanying manuscript (including the optional *inter-annual state noise* used to reproduce the “noisy” simulation scenarios).

The main files to run are pertexp_GUI.app, pertexp_simulator.m and RC_contour.
---

## Repository contents

| File | Type | Purpose |
|---|---|---|
| `pertexp_simulator.m` | **Script** | Standalone simulator **without GUI**. Runs a 5-year simulation (weekly steps), optionally adds inter-annual state noise, computes the age-specific annual cumulative infections and vaccinations, and displays and saves the corresponding figures. |
| `pertexp_GUI.mlapp` | **MATLAB App** | Graphical user interface to set vaccination parameters and initial conditions, optionally enable noise, and display interactive bar plots (annual, age-specific) over a **fixed 5-year horizon**. |
| `pertexp_run.m` | **Function** | Core routine called by the GUI: runs the simulation and returns the annual quantities displayed by the interactive bar plots. |
| `discrete_system.m` | **Function** | The **discrete-time model** (weekly step). Computes the weekly state trajectory and weekly flows (incidence, vaccinations, deaths, etc.). |
| `add_noise_state.m` | **Function** | Adds **multiplicative Gaussian noise** to the state vector at yearly boundaries (inter-annual variability), with renormalisation by age class. |
| `RC_calc.m` | **Function** | Computes the **disease-free equilibrium (DFE)** and the reproduction numbers \(\mathcal{R}_0\) and \(\mathcal{R}_c\). |
| `RC_contour.m` | **Script** | Generates contour plots of \(\mathcal{R}_c\) as a function of vaccination coverages (maternal vs infant or maternal vs booster). |

---

## Requirements

- MATLAB with App Designer support (to run `pertexp_GUI.mlapp`).

No additional toolboxes are required beyond standard MATLAB functionality.

---

## Quick start

1. Clone/download the repository.
2. Open MATLAB and set the **Current Folder** to the repository root.

---

## Citation

If you use these codes in academic work, please cite the corresponding PerTexP manuscript.
