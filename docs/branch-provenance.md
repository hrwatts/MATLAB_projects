# Branch Provenance

This branch reorganizes branch snapshots from the source repository [`hrwatts/MATLAB_projects`](https://github.com/hrwatts/MATLAB_projects) into a single folder-based layout in the fork `hrwdata/MATLAB_projects`.

The imports are snapshot-based and intentionally self-contained. Repeated files were left inside the project folders where they originally appeared.

| Source branch | Source commit | Destination folder | Notes |
| --- | --- | --- | --- |
| `main` | `f46c5e3f9920a11d215464230296339db5fc9c43` | `projects/eulers-method/` | Imported `eulers.m` only; root README replaced by repository index |
| `tutorial` | `e4509e93fbd25d08e43c4bdd260e2823353a4cdb` | `projects/tutorial/` | `.html` and `.pdf` exports moved into `exports/` |
| `bifurcation` | `4a664825be4816ead316e0cb3894a27bfddcf0f8` | `projects/bifurcation/` | Original `README-b.md` preserved as `README.source.md` |
| `SIR-COVID` | `c1e85e01a54dd8a830418a9fc8da94f3bf2126ca` | `projects/sir-covid/` | Original `README-SIR.md` preserved as `README.source.md` |
| `gradient_descent` | `4170861e4e8b8e3cc0a2893893f13546635f06bf` | `projects/gradient-descent/` | Preserves local copy of `eulers.m` |
| `simplex_method` | `8e77cdb2c1366886a2345880681ba48cb6b5b77c` | `projects/simplex-method/` | Preserves local copy of `eulers.m` |

## Import Policy

- Each branch was imported into its own folder as a snapshot, not merged into the repository root.
- Original branch readmes were preserved as `README.source.md` where available.
- Large tutorial exports were preserved and separated into `projects/tutorial/exports/`.
