# SIR-COVID

This project fits a simple SIR model to COVID case data and includes the data-cleaning and weekly aggregation steps used to support that fit.

## Files

- `SIR_fit_AL.mlx`: live script workflow for the Alabama dataset
- `al_covid.csv`: case and death data used by the live script
- `week_tab.m`: converts raw data into weekly aggregates
- `sir_fit.m`: fits an SIR curve and estimates parameters
- `sir_model.m`: ODE right-hand side for the SIR system
- `newton_m.m`: Newton iteration helper used during parameter estimation
- `README.source.md`: original branch README

## Usage

Open `SIR_fit_AL.mlx` in MATLAB for the end-to-end workflow, or call `sir_fit` directly with a CSV file path and population value:

```matlab
[tspan, curve, k, a] = sir_fit("al_covid.csv", 227900);
```
