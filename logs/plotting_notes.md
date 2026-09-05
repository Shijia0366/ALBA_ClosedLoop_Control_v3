# Plotting warning record

The final `v3_post.py` plotting/export run completed with exit code 0.
Matplotlib emitted this non-fatal performance warning at the poster export:

```text
UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
```

This concerns automatic legend placement speed, not the simulation solver or data.
PNG, SVG and PDF files were generated. The final poster PNG was visually inspected.
MATLAB/Simulink warnings and execution metadata remain in each original case MAT,
the case run logs, and `environment_and_build.log` / `final_wiring_rerun.log`.
