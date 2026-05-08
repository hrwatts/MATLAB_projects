# Algorithms

## Numerical ODE Methods

- `src/numerics/eulers.m`: Euler method for first-order initial value problems
- `src/dynamics/plotd.m`: normalized direction-field plotting helper
- `src/dynamics/plotb.m`: bifurcation-style plotting helper built on repeated Euler simulations

## Optimization

### Line Search

- `src/optimization/line_search/armijo.m`
- `src/optimization/line_search/goldstein.m`
- `src/optimization/line_search/wolfe.m`

These functions implement line search rules for optimization routines that work with symbolic objective functions.

### Conjugate Gradient

- `src/optimization/conjugate_gradient/pcgrad.m`

This routine implements a partial conjugate gradient method with restart behavior controlled by `m`.

### Linear Programming

- `src/optimization/linear_programming/bounded_simplex.m`
- `src/optimization/linear_programming/primal_dual.m`

These routines contain prototype research implementations of bounded simplex and primal-dual methods.

## Epidemiology

- `src/epidemiology/week_tab.m`: aggregates cumulative case/death data into weekly totals
- `src/epidemiology/sir_model.m`: SIR right-hand side
- `src/epidemiology/sir_fit.m`: parameter-fitting workflow using weekly data and ODE simulation
