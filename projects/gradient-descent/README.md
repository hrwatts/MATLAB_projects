# Gradient Descent

This folder contains optimization routines from the original `gradient_descent` branch, centered on a partial conjugate gradient implementation and several line-search strategies.

## Files

- `pcgrad.m`: partial conjugate gradient method
- `armijo.m`: Armijo line search
- `wolfe.m`: Wolfe line search
- `goldstein.m`: Goldstein line search
- `eulers.m`: branch-local helper copy preserved to keep the imported snapshot self-contained
- `README.source.md`: original branch README

## Usage

`pcgrad.m` expects a symbolic scalar objective and an initial point. A minimal call shape is:

```matlab
syms x1 x2
f = (x1 - 1)^2 + (x2 + 2)^2;
[x_star, z_star, x_values, iterations] = pcgrad(f, [0; 0], 2, 'A', 'F', 1e-12, 50);
```

Use `m = 0` for steepest descent behavior and `m = n` for conjugate-gradient-style restarts, matching the header comments in `pcgrad.m`.
