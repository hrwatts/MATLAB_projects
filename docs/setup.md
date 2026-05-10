# Setup

## MATLAB Requirements

- Base MATLAB
- ODE solvers for epidemiology examples
- Symbolic Math Toolbox for:
  - `src/optimization/conjugate_gradient/pcgrad.m`
  - `src/epidemiology/sir_fit.m`

## Recommended Startup

From the repository root:

```matlab
addpath(genpath('src'));
addpath(genpath('examples'));
addpath(genpath('sample-data'));
```

## Layout Notes

- Use `src/` for reusable functions only.
- Use `examples/` for runnable demos and Live Scripts.
- Use `tests/` for `matlab.unittest` coverage.
- Keep the preserved `projects/` tree only as a temporary migration aid during this first reorganization pass.
