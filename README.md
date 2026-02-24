# Heterogeneous Transfer Learning for High-dimensional Regression with Feature Mismatch
This repository provides simulation and data example codes for the paper https://arxiv.org/abs/2412.18081

## Structure

- **Simulations/** - Comprehensive simulation comparing HTL, HmTL, and baseline methods
  - See [Simulations/README.md](Simulations/README.md) for setup and usage
  - Quick start: [Simulations/QUICKSTART.md](Simulations/QUICKSTART.md)
  
- **CuratedOvarianCancer/** - Real data application example
  - See [CuratedOvarianCancer/REPRODUCIBILITY.md](CuratedOvarianCancer/REPRODUCIBILITY.md)

## Key Files

- `sv-cv.R` - Sieve estimation with cross-validation (dependency for simulations)
- `sv-prim.R` - Utility functions for non-linear features (dependency for simulations)
