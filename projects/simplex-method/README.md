# Simplex Method

This project contains bounded simplex and primal-dual linear optimization routines together with example scripts.

## Files

- `bounded_simplex.m`: bounded simplex solver
- `bounded_simplex_example.m`: sample bounded simplex problems
- `primal_dual.m`: primal-dual solver
- `primal_dual_example.m`: sample primal-dual problems
- `eulers.m`: branch-local helper copy preserved to keep the imported snapshot self-contained
- `README.source.md`: original branch README

## Usage

Run the example scripts from this folder in MATLAB:

```matlab
bounded_simplex_example
primal_dual_example
```

Both examples define matrices and vectors inline, call the corresponding solver, and pause between example cases.
