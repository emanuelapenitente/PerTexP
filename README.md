# PerTexP (Pertussis Time Exploration)
PerTexP is a MATLAB-based, interactive modelling tool designed to investigate the pertussis transmission dynamics and to support the evaluation of different vaccination strategies and short-term projections. 

The implementation of PerTexP is built upon a discrete-time compartmental model with the two age classes, *infants*, comprising individuals from 0 to 11 months, and *non-infants*, comprising individuals aged 1 year and older. All the details on the model formulation are reported in the accompanying manuscript.

This repository contains the MATLAB scripts and functions to **simulate different vaccination scenarios** and to compute the control reproduction number $\mathcal{R}_c$ under different vaccination strategies. This repository has three main files:
  
- `PerTexP_GUI.mlapp` — a MATLAB app with a Graphical User Interface (GUI), built with MATLAB App Designer, which allows users to interactively explore and evaluate different vaccination strategies in the context of the 2024 Italian pertussis outbreak,

- `pertexp_simulator.m` — run simulations without the GUI,

- `RC_contour.m` — generate contour plots of the control reproduction number $\mathcal{R}_c$ as a function of vaccination coverages.

**Note:** PerTexP is **applied to the 2024 Italian pertussis outbreak**, which is used in the manuscript as an *illustrative case study* (demography/initial conditions and baseline parameter values reflect the Italian context). This configuration is **not meant to be restrictive**: PerTexP is a general modelling framework and can be adapted to **any country/setting** by updating the demographic/epidemiological inputs and initial conditions.


## Repository contents

| File | Type | Purpose |
|---|---|---|
| `PerTexP_GUI.mlapp` | **MATLAB App** | Graphical User Interface to set vaccination parameters and initial conditions, optionally enable noise, and display interactive bar plots (annual, age-specific) over a **fixed 5-year horizon**. |
| `pertexp_logo.png` | **Image asset** | PerTexP logo used by the GUI (`pertexp_GUI.mlapp`) for display purposes (has a merely aesthetic purpose, it is not required to run the simulator/analysis scripts). |
| `pertexp_run.m` | **Function** | Core routine called by the GUI: runs the simulation and returns the annual quantities displayed by the interactive bar plots. |
| `pertexp_simulator.m` | **Script** | Standalone simulator **without GUI**. Runs a 5-year simulation (weekly steps), optionally adds inter-annual state noise, computes the age-specific annual cumulative infections and vaccinations, and displays and saves the corresponding figures. |
| `discrete_system.m` | **Function** | The **discrete-time model** (weekly step). Computes the weekly state trajectory and weekly flows (incidence, vaccinations, deaths, etc.). |
| `add_noise_state.m` | **Function** | Adds **multiplicative Gaussian noise** to the state vector at yearly boundaries (inter-annual variability), with renormalisation by age class. |
| `RC_calc.m` | **Function** | Computes the **disease-free equilibrium (DFE)** and the reproduction numbers $\mathcal{R}_0$ and $\mathcal{R}_c$. |
| `RC_contour.m` | **Script** | Generates contour plots of $\mathcal{R}_c$ as a function of vaccination coverages (maternal vs infant or maternal vs booster). |

## Requirements

- MATLAB with App Designer support (to run `pertexp_GUI.mlapp`),
- The scripts use `exportgraphics` to save figures (built-in, no toolbox), which requires MATLAB R2020a or newer. If using older releases, replace `exportgraphics` calls with `print` or `saveas`.

No additional toolboxes are required.

## Usage

> **Important:** keep all repository files in the **same folder** (i.e., do not move `.m` files, the `.mlapp`, or `pertexp_logo.png` into subfolders). The GUI and the scripts assume that all contents are located in a single directory; moving files may break function calls.

### 1) How to use the GUI

Open and run `PerTexP_GUI.mlapp`.

#### What you can change in the app

- Vaccination parameters (vaccination coverages, efficacies, upatakes and duration of immunity),
- Initial conditions of the simulation,
- Optional inter-annual state noise.

#### What the app displays

The GUI displays the **value of the control reproduction number**, $\mathcal{R}_c$, and **interactive bar plots** over a **fixed 5-year horizon**, including:

- Annual cumulative **incidence** (age-specific),
- Annual cumulative **vaccinations/boosters** (age-specific). Infant vaccinations are reported as a percentage of the total infant population; non-infant cumulative annual boosters are reported over 10,000 of the total non-infant population (both infant and non-infant populations are kept fixed in the time-horizon). 

Internally, the app calls:
- `pertexp_run.m` (runs the model and computes annual quantities),
- `pertexp_run.m` in turn calls `discrete_system.m` (weekly simulation), `RC_calc.m` (computes the control reproduction number) and, if enabled, `add_noise_state.m` at year boundaries.

### 2) How to use the PerTexP simulator (without the GUI)

Run the main script ```pertexp_simulator.m```.

#### What it does

- Sets parameters and initial conditions,
- Simulates the model for **5 years** with a **weekly time step**,
- Computes annual relevant quantities (age-specific cumulative annual incidence and vaccinations),
- Produces and saves the figures in both .png and .eps formats.

#### Noise option (inter-annual variability)

In `pertexp_simulator.m`, the *Simulation* section contains:
```matlab
use_noise = false;   % true = add noise, false = no noise

if use_noise
    noise_level = 0.05;   % ±5%
else
    noise_level = 0;
end
```

- If `use_noise = true`, the script perturbs the **end-of-year state** before starting the next year using:
  - `add_noise_state(x_end, N1, N2, noise_level)`.
- The default noise magnitude is `noise_level = 0.05` (i.e., ~5% multiplicative perturbations, clamped to avoid extreme values).
To obtain reproducible noisy runs, set a random seed before running the script, e.g.:

```matlab
rng(1);
pertexp_simulator
```

#### Outputs

`pertexp_simulator.m` creates (if missing) and writes figures to:

- `fig_eps/` (EPS),
- `fig_png/` (PNG).

### 3) How to produce the contour plots of $\mathcal{R}_c$

Run ```RC_contour.m```.

This script produces contour plots of $\mathcal{R}_c$ over grids of vaccination coverages, including:

1. $\mathcal{R}_c$ as a function of **annual infant coverage** $\psi_1^{\mathrm{ann}}$ and **maternal coverage** $p$,
2. $\mathcal{R}_c$ as a function of **annual booster coverage** $\psi_2^{\mathrm{ann}}$ and **maternal coverage** $p$.

Figures are saved to:
- `fig_contour_eps/` (EPS),
- `fig_contour_png/` (PNG).

## Citation

If you use these codes in academic work, please cite the corresponding PerTexP manuscript.

## Contact

Suggestions and feedback are very welcome — feel free to open an issue on GitHub or get in touch.
