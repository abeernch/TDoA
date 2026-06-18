# TDoA / FDoA

MATLAB research and test scripts for time-difference of arrival (TDoA),
frequency-difference of arrival (FDoA), hyperbola and hyperboloid construction,
geometry analysis, and related geolocation experiments.

The repository is organized around analysis workflows rather than a packaged
application. It contains:

- 2D and 3D TDoA solver and plotting experiments.
- GDoP and CRLB-related geometry analysis.
- Least-squares and benchmarking scripts.
- Static-fusion and integration test harnesses.
- Shared utility functions reused across script families.

## Repository Layout

| Folder | Purpose |
| --- | --- |
| `220601_Static Fusion TDOA/` | Static fusion experiments and supporting scripts. |
| `220901_2D Custom Hyperbolae Plotting and Solution/` | Hyperbola construction, coordinate transforms, and 2D plotting experiments. |
| `221201_ELINTSG TDoA Integration/` | Integration scripts for TOA/TDoA experiments and 3D adaptation. |
| `250103_GDOP/` | GDOP and covariance-related functions and demos. |
| `250110_TDoA_Benchmarking/` | Benchmarking scripts and helper functions. |
| `250129_LeastSquares_Custom_TDoA/` | Least-squares TDoA scripts and supporting functions. |
| `Functions/` | Shared utility functions used by multiple script sets. |
| `Test Scripts/` | General verification scripts for GDoP, CRLB, and related behaviors. |
| `Test Scripts TDOA/` | TDoA-specific test scripts and auxiliary examples. |

## Public Documentation Scope

This repository keeps the public documentation generic. Deployment-specific site
names, coordinates, and similar location identifiers are intentionally omitted
from the root documentation.

## Getting Started

1. Open the repository in MATLAB.
2. Add the repository folder and subfolders to the MATLAB path.
3. Open the script family you want to inspect or run.
4. Start with the test scripts if you want to verify the math before working
   through the benchmark or integration folders.

## Representative Entry Points

- `test_gridSearchSolplot_241230.m`
- `test_fdoa_isocontours_241007.m`
- `test_fdoa_MLE_250127.m`
- `test_coarsegridSearchSolplot_250101.m`
- `250110_TDoA_Benchmarking/tdoaGDoP_bench_250110.m`
- `250129_LeastSquares_Custom_TDoA/mainadsb.m`

## Related Paper

For the complementary geometry-quality metric work that matches this codebase,
cite:

```text
Abeer Nasir Chaudhry, Salman Liaquat, and Muhammad Mohsin Khadim,
$\kappa$: A Geometry-Quality Metric Complementary to GDoP for Closed-Form
TDoA Multilateration, arXiv:2606.14372, 2026.
```

## Notes

- The repository contains MATLAB `.m` scripts and a few live scripts (`.mlx`).
- Several directories contain archived iterations of the same analysis idea.
- Some helper files are shared across multiple experiment folders.

## License

Released under the MIT License. See [LICENSE](LICENSE).
