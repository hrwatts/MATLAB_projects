# Euler's Method

This folder contains the base Euler method routine from the original `main` branch.

## Files

- `eulers.m`: computes Euler approximations on an interval `[a,b]` and returns the sampled `x` values, approximate solution values, slopes, and a MATLAB table

## Usage

Call the function with a derivative handle and initial value problem data:

```matlab
[x, y, f, table] = eulers(@(x,y) x + y, 0, 1, 1, 0.1);
```

This implementation is also reused in several of the other imported projects, but it is kept here as a standalone showcase of the original `main` branch content.
