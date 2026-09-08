# ALBA Step 9 — Closed-Loop Control Simulation v3

MATLAB/Simulink implementation and numerical validation of the closed-loop
control study developed for the ALBA soft robotic bladder-assist system.

This repository supports the modelling and control results reported in the MSc
thesis:

**Development of an Integrated Sensing and Control Platform for a Soft Robotic
Bladder Assist System**

## Scope

The repository contains:

- the nonlinear motor-driven ALBA reference model;
- MATLAB/Simulink closed-loop implementation;
- fixed-rate, pure-PI and geometric feedforward-PI controller comparisons;
- startup and volume-based termination studies;
- frozen-controller simulations for 200, 300 and 400 mL targets;
- reduced-order/full-model comparisons;
- independent Python-to-Simulink numerical verification;
- simulation results and validation records.

The closed-loop results are simulation-based. Physical experiments reported in
the thesis were performed separately for fixed-load discharge characterisation.

## Final configuration

The final controller was selected using the 300 mL case and then frozen for
the 200, 300 and 400 mL evaluations.

| Parameter | Final value |
|---|---:|
| Flow reference | 15 mL/s |
| Controller | Pure PI |
| `Kp` | 50 |
| `Ki` | 600 |
| `Kaw` | 10 |
| Startup duration | 3 s |
| Terminal stop threshold | 0.05 mL |
| Command limit | 0–3000 STEP/s |
| Nominal hydrostatic head | +0.10 m |

Pure PI was retained because it satisfied the plateau-flow requirement while
giving a lower simulated peak torque than geometric feedforward-PI in the
final controller comparison.

## Final simulations

| Target volume | Plateau flow RMSE | Maximum relative deviation | Completion time |
|---:|---:|---:|---:|
| 200 mL | 0.0304 mL/s | 0.494% | 18.55 s |
| 300 mL | 0.0393 mL/s | 0.494% | 25.20 s |
| 400 mL | 0.0619 mL/s | 0.941% | 31.83 s |

All three nominal simulations satisfied the predefined ±5% plateau-flow
criterion.

For the 300 mL case, the simulated peak torque was approximately 0.2808 N·m,
below the conditional model torque envelope of approximately 0.2966 N·m.

These results apply to the nominal model assumptions and ideal simulated
feedback; they do not constitute hardware closed-loop validation.

## Repository structure

```text
ALBA_Step9_ClosedLoop_v3/
│
├── README.md
├── requirements-v3.txt
├── ALBA_Stage8_Results.mat
│
├── rebuild_v3.m
├── run_v3_case.m
├── run_v3_final_validation.m
├── run_v3_optimisation_screen.m
├── alba9_*.m
├── v3_*.m
├── v3_*.py
│
├── stage9b_port/
│   └── alba9b_full_loop.slx
│
├── frozen_reference/
├── frozen_rom/
├── frozen_v2_ROM_results/
├── provenance_v2/
│
├── results/
│   ├── D_diagnostic_baseline/
│   ├── D_convergence_tight/
│   ├── A_startup_2p5/
│   ├── A_startup_3/
│   ├── B_terminal_0p05/
│   ├── B_terminal_0p25/
│   ├── B_guard_localised/
│   ├── C_combined_300/
│   ├── FINAL_verified/
│   └── FINAL_verified_tight/
│
├── validation/
└── logs/
```

The principal nonlinear closed-loop Simulink model is:

```text
stage9b_port/alba9b_full_loop.slx
```

The final frozen-controller results reported in the thesis are stored in:

```text
results/FINAL_verified/
```

## Reproduction

The reported simulations were executed using MATLAB/Simulink R2025b.

Install the Python dependencies:

```bash
python -m pip install -r requirements-v3.txt
```

Then run from the repository root in MATLAB:

```matlab
rebuild_v3('C:\path\to\python.exe')
```

`rebuild_v3.m` rebuilds and runs the Simulink cases, performs the numerical
verification, and regenerates the derived validation outputs.

Run the rebuild in a separate copy of the repository if the archived result
files are to remain unchanged.
