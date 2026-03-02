# PerTexP (Pertussis Time Exploration)
PerTexP is a MATLAB-based, interactive modelling tool designed to investigate pertussis transmission dynamics and to support the evaluation of vaccination strategies and short-term projections. 
This repository contains the MATLAB scripts and functions to **simulate different vaccination scenarios** and to **compute the control reproduction number** $\mathcal{R}_c$ under different vaccination strategies.

The implementation of PerTexP is based on an underlying **discrete-time epidemic model** described in the accompanying manuscript (including the optional *inter-annual state noise* used to reproduce the “noisy” simulation scenarios).

## Repository contents

| File | Type | Purpose |
|---|---|---|
| `pertexp_simulator.m` | **Script** | Standalone simulator **without GUI** (Graphical User Interface). Runs a 5-year simulation (weekly steps), optionally adds inter-annual state noise, computes the age-specific annual cumulative infections and vaccinations, and displays and saves the corresponding figures. |
| `pertexp_GUI.mlapp` | **MATLAB App** | Graphical user interface to set vaccination parameters and initial conditions, optionally enable noise, and display interactive bar plots (annual, age-specific) over a **fixed 5-year horizon**. |
| `pertexp_run.m` | **Function** | Core routine called by the GUI: runs the simulation and returns the annual quantities displayed by the interactive bar plots. |
| `discrete_system.m` | **Function** | The **discrete-time model** (weekly step). Computes the weekly state trajectory and weekly flows (incidence, vaccinations, deaths, etc.). |
| `add_noise_state.m` | **Function** | Adds **multiplicative Gaussian noise** to the state vector at yearly boundaries (inter-annual variability), with renormalisation by age class. |
| `RC_calc.m` | **Function** | Computes the **disease-free equilibrium (DFE)** and the reproduction numbers $\mathcal{R}_0$ and $\mathcal{R}_c$. |
| `RC_contour.m` | **Script** | Generates contour plots of $\mathcal{R}_c$ as a function of vaccination coverages (maternal vs infant or maternal vs booster). |

## Requirements

- MATLAB with App Designer support (to run `pertexp_GUI.mlapp`).
- The scripts use `exportgraphics` to save figures (built-in, no toolbox), which requires MATLAB R2020a or newer.  
  If using older releases, replace `exportgraphics` calls with `print` or `saveas`.

No additional toolboxes are required.

## Usage

This repository has three main files:

- `pertexp_simulator.m` — run simulations without the Graphical User Interface (GUI)

- `RC_contour.m` — generate contour plots of the control reproduction number $\mathcal{R}_c$
- `pertexp_GUI.mlapp` — interactively explore scenarios via a GUI (App Designer)

### 1) How to use the GUI (Graphical User Interface)

Open and run `pertexp_GUI.mlapp`.

#### What you can change in the app

- Vaccination parameters (vaccination coverages, efficacies, upatakes and duration of immunity).
- Initial conditions of the simulation.
- Optional inter-annual state noise.

#### What the app displays

The GUI displays the **value of the control reproduction number**, $\mathcal{R}_c$, and **interactive bar plots** over a **fixed 5-year horizon**, including:

- Annual cumulative **incidence** (age-specific).
- Annual cumulative **vaccinations/boosters** (age-specific). Infant vaccinations are reported as a percentage of the total infant population; non-infant cumulative annual boosters are reported over 10,000 of the total non-infant population (both infant and non-infant populations are kept fixed in the time-horizon). 

Internally, the app calls:
- `pertexp_run.m` (runs the model and computes annual quantities)
- `pertexp_run.m` in turn calls `discrete_system.m` (weekly simulation), `RC_calc.m` (computes the control reproduction number) and, if enabled, `add_noise_state.m` at year boundaries.

### 2) How to use the PerTexP simulator (without the GUI)

Run the main script ```pertexp_simulator.m```.

#### What it does

- Sets parameters and initial conditions.
- Simulates the model for **5 years** with a **weekly time step**.
- Computes annual relevant quantities (age-specific cumulative annual incidence and vaccinations).
- Produces and saves the figures in both .png and .eps formats.

#### Noise option (inter-annual variability)

In `pertexp_simulator.m`, the *Simulation* section contains:
```matlab
use_noise = false;   % true = add noise, false = no noise
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


Run: ```RC_contour.m```

This script produces contour plots of \(\mathcal{R}_c\) over grids of vaccination coverages, including:

1. $\mathcal{R}_c$ as a function of **annual infant coverage** $\psi_1^{\mathrm{ann}}$ and **maternal coverage** $p$.
2. $\mathcal{R}_c$ as a function of **annual booster coverage** $\psi_2^{\mathrm{ann}}$ and **maternal coverage** $p$.

Figures are saved to:
- `fig_contour_eps/` (EPS),
- `fig_contour_png/` (PNG).

## Citation

If you use these codes in academic work, please cite the corresponding PerTexP manuscript.

## Contact

Suggestions and feedback are very welcome — feel free to open an issue on GitHub or get in touch.
