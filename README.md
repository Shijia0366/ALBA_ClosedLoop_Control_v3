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
- startup-time evaluation;
- volume-based termination analysis;
- frozen-controller simulations for 200, 300 and 400 mL targets;
- reduced-order/full-model comparisons;
- independent Python-to-Simulink numerical verification;
- simulation results and reproducibility records.

The closed-loop results in this repository are simulation-based. Physical
experiments in the thesis were performed separately for fixed-load discharge
characterisation.

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
| Terminal stop threshold | 0.05 mL before closure threshold |
| Command limit | 0–3000 STEP/s |
| Nominal hydrostatic head | +0.10 m |

Pure PI was retained because it satisfied the flow-tracking requirement while
giving a lower simulated peak torque than the geometric feedforward-PI
controller in the final controller comparison.

## Final simulations

The frozen controller was evaluated at three target discharged volumes.

| Target volume | Plateau flow RMSE | Maximum relative deviation | Completion time |
|---:|---:|---:|---:|
| 200 mL | 0.0304 mL/s | 0.494% | 18.55 s |
| 300 mL | 0.0393 mL/s | 0.494% | 25.20 s |
| 400 mL | 0.0619 mL/s | 0.941% | 31.83 s |

All three nominal simulations satisfied the predefined ±5% plateau-flow
criterion.

For the 300 mL case, the simulated peak torque was approximately
0.2808 N·m, below the conditional model torque envelope of approximately
0.2966 N·m.

## Repository structure

```text
ALBA_Step9_ClosedLoop_v3/
│
├── README.md
├── REPORT_v3.md
├── CHANGELOG_v3.md
├── MANIFEST_v3.json
├── requirements-v3.txt
│
├── rebuild_v3.m
├── run_v3_case.m
├── alba9_*.m
│
├── stage9b_port/
│   ├── alba9b_full_loop.slx
│   └── ...
│
├── frozen_reference/
├── frozen_rom/
│
├── results/
│   ├── D_diagnostic_baseline/
│   ├── A_startup_*/
│   ├── B_terminal_*/
│   ├── B_guard_localised/
│   ├── C_combined_300/
│   ├── FINAL_verified/
│   └── FINAL_verified_tight/
│
├── validation/
│   ├── python_trajectories/
│   ├── ROM_full_*/
│   └── ...
│
├── poster/
└── baseline_v2_reports/
