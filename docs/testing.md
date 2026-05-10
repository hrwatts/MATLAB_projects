# Testing

## Test Strategy

The `tests/` directory contains `matlab.unittest` classes for:

- numerical utilities
- epidemiology helpers
- optimization routines
- smoke tests for example-facing behavior

## Running Tests

```matlab
addpath(genpath('src'));
results = runtests('tests');
table(results)
```

## Toolbox Notes

- Tests for `pcgrad` and `sir_fit` require Symbolic Math Toolbox.
- If Symbolic Math Toolbox is unavailable, those tests should be skipped rather than failing the whole suite.

## Current Limitations

- MATLAB execution could not be validated in the current environment during this reorganization pass because the local MATLAB license is not usable.
