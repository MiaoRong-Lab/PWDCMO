# PlatEMO Adapter

This folder contains files arranged in the same style as the official PlatEMO project.

To use them without copying files, add both official PlatEMO and this folder to the MATLAB path:

```matlab
addpath(genpath('path/to/PlatEMO/PlatEMO'));
addpath(genpath('path/to/PWDCMO/PlatEMO'));
```

For an upstream PlatEMO contribution, copy or submit:

- `Algorithms/Multi-objective optimization/PWDCMO/PWDCMO.m`
- `Problems/Multi-objective optimization/PWDCMO-DCMOP/PWDCMO_DCMOP.m`
- `Problems/Multi-objective optimization/PWDCMO-DCMOP/PWDCMO_DCMOPEvaluate.m`
- `Problems/Multi-objective optimization/PWDCMO-DCMOP/PWDCMO_DCMOPPF.m`

The problem class is parameterized by `instance` from 1 to 8. If preferred by the PlatEMO maintainers, it can be split into eight separate problem classes later.
