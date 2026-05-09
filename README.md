# PWDCMO

[![DOI](https://img.shields.io/badge/DOI-10.1109%2FTEVC.2024.3418470-blue.svg)](https://doi.org/10.1109/TEVC.2024.3418470)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2018a%2B-orange.svg)](#requirements)
[![PlatEMO](https://img.shields.io/badge/PlatEMO-compatible-green.svg)](#platemo-usage)

MATLAB code for the paper:

> D. Gong, M. Rong, N. Hu, Y. Wang, W. Pedrycz, and S. Yang, "A Prediction and Weak Coevolution-Based Dynamic Constrained Multiobjective Optimization," IEEE Transactions on Evolutionary Computation, 2025, 29(4): 1328-1342. DOI: [10.1109/TEVC.2024.3418470](https://doi.org/10.1109/TEVC.2024.3418470).

Repository: <https://github.com/MiaoRong-Lab/PWDCMO>

This repository contains two entry points:

- `Standalone/`: standalone reproduction code derived from the original experiment scripts.
- `scripts/run_paper_experiment.m`: standalone reproduction script based on the original experiment code.
- `PlatEMO/`: a PlatEMO-style adapter containing the `PWDCMO` algorithm and `PWDCMO_DCMOP` dynamic constrained test problem.

## Requirements

- MATLAB R2018a or later.
- For the PlatEMO adapter: PlatEMO v4.x or later.

The standalone path no longer requires `normrnd` or `pdist2`; the refactor replaced those calls with base MATLAB equivalents. The hypervolume helper still contains a Monte Carlo branch for four or more objectives, but the provided DCMOP instances are bi-objective.

## Repository Layout

```text
PWDCMO/
  Standalone/   Original-style code for reproducing the paper experiments
  PlatEMO/      Files that can be copied into or added alongside PlatEMO
  scripts/      Small entry scripts for smoke tests and experiments
  README.md
  CITATION.cff
```

## Standalone Usage

Run a quick check:

```matlab
run('scripts/run_smoke_test.m')
```

Run the full standalone experiment:

```matlab
run('scripts/run_paper_experiment.m')
```

The full script defaults to DCMOP1, `N = 200`, `pps = 100`, 20 independent runs, and time points `0:0.25:20`. Results are saved under `results/`.

## PlatEMO Usage

Copy or add the repository's `PlatEMO` folder to the MATLAB path together with PlatEMO:

```matlab
addpath(genpath('path/to/PlatEMO/PlatEMO'));
addpath(genpath('path/to/PWDCMO/PlatEMO'));

Problem   = PWDCMO_DCMOP('N',100,'D',10,'maxFE',20000,'parameter',{1,20,4});
Algorithm = PWDCMO('save',1,'metName',{'IGD','HV'});
Algorithm.Solve(Problem);
```

`PWDCMO_DCMOP` parameters are:

- `instance`: DCMOP instance number from 1 to 8.
- `taut`: number of generations for one environment.
- `nt`: number of time steps per unit time. The default `nt = 4` gives time increments of `0.25`.

## Files

- `Standalone/main.m`: standalone experiment loop.
- `Standalone/NSGA2_1.m`: static optimizer with constrained and diversity archives.
- `Standalone/dynamic_response_1.m`: prediction-based response after an environmental change.
- `Standalone/variable_move.m`: infeasible solution movement toward feasible nondominated guidance solutions.
- `PlatEMO/Algorithms/Multi-objective optimization/PWDCMO/PWDCMO.m`: PlatEMO-compatible algorithm class.
- `PlatEMO/Problems/Multi-objective optimization/PWDCMO-DCMOP/PWDCMO_DCMOP.m`: PlatEMO-compatible dynamic constrained problem class.

## License

This repository is released under the research-use license in `LICENSE`,
following the research and educational use spirit of PlatEMO. Commercial use
requires separate permission from the copyright holders.

## Notes Before Public Release

- Validate the PlatEMO adapter against the final paper tables before submitting it upstream.
- If this code or the PlatEMO adapter is used with PlatEMO, acknowledge and cite PlatEMO as required by the PlatEMO project.

## Citation

```bibtex
@article{Gong2025PWDCMO,
  title   = {A Prediction and Weak Coevolution-Based Dynamic Constrained Multiobjective Optimization},
  author  = {Gong, Dunwei and Rong, Miao and Hu, Na and Wang, Yan and Pedrycz, Witold and Yang, Shengxiang},
  journal = {IEEE Transactions on Evolutionary Computation},
  volume  = {29},
  number  = {4},
  pages   = {1328--1342},
  year    = {2025},
  doi     = {10.1109/TEVC.2024.3418470}
}
```

When using PlatEMO, also cite:

```bibtex
@article{PlatEMO,
  title   = {{PlatEMO}: A {MATLAB} platform for evolutionary multi-objective optimization},
  author  = {Tian, Ye and Cheng, Ran and Zhang, Xingyi and Jin, Yaochu},
  journal = {IEEE Computational Intelligence Magazine},
  volume  = {12},
  number  = {4},
  pages   = {73--87},
  year    = {2017}
}
```
