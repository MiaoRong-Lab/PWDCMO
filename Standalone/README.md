# Standalone Code

This folder contains the original-style MATLAB implementation for reproducing the PWDCMO experiments without installing PlatEMO.

Use the scripts in the repository root:

```matlab
run('scripts/run_smoke_test.m')
run('scripts/run_paper_experiment.m')
```

`main.m` is the standalone experiment loop. The script `scripts/run_paper_experiment.m` sets the paper-style defaults and saves results under `results/`.
