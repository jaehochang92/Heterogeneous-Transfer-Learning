# Heterogeneous Transfer Learning Simulation

## Overview

This simulation compares four machine learning approaches for transfer learning scenarios where source and target data may have different feature spaces:

1. **HTL-Linear (Heterogeneous Transfer Learning with Linear Mapping)**
   - Learns linear mappings from source features to missing target features
   - Uses the linear feature mappings combined with target features in a lasso model
   
2. **HTL-Sieve (Heterogeneous Transfer Learning with Sieve Estimation)**
   - Learns non-linear feature mappings via sieve basis estimation
   - Uses complex basis functions to capture non-linearities in feature relationships
   - Combined with target features in a lasso model

3. **HmTL (Homogeneous Transfer Learning)**
   - Uses source data that shares the same features as target data
   - Combines source and target data directly for training
   - Baseline for scenarios with feature alignment across domains

4. **Lasso-Target-Only**
   - Standard lasso regression using only target domain training data
   - Primary baseline for comparison

## Files Description

### Core Simulation Files (in Simulations/ directory)

- **`simulation.R`**: Main simulation code containing:
  - Data generation functions for heterogeneous/homogeneous scenarios
  - Implementation of all four methods
  - Main simulation loop that runs multiple iterations
  - Performance metrics (RMSE, MAE)

- **`run_simulation.R`**: Convenient wrapper to run simulations with different configurations
  - Loads all dependencies and provides a `run_simulation()` function
  - Examples for small, medium, and large simulations

- **`summarize_sim.R`**: Statistical summary and comparison script
  - Outputs summary statistics by method
  - Calculates improvement percentages over baseline
  - Performs pairwise t-tests
  - Can be run from command line: `Rscript summarize_sim.R results.rds`

- **`plot_sim_results.R`**: Visualization script
  - Creates boxplots, violin plots, and bar charts
  - Saves publication-ready plots
  - Can be run from command line: `Rscript plot_sim_results.R results.rds output.pdf`

### Dependencies (in parent directory)

- **`sv-cv.R`** (parent dir): Sieve estimation with cross-validation (from Zhang & Simon, 2023)
  - `cross_validated_sieve()`: CV for sieve basis size and regularization
  - `fit_sieve_with_cv()`: Fit sieve model with optimal hyperparameters
  - `cv_fitH_from_prxy()`: Fit feature mappings (linear or sieve) in parallel
  - `fmap_predict()`: Predict using learned feature mappings

- **`sv-prim.R`** (parent dir): Utility functions for non-linear feature transformations
  - `H()`, `h()`, `hj()`: Additive non-linear basis functions for testing

## Simulation Parameters

When running `main()` or `run_simulation()`, you can customize:

```r
main(
  n_sim = 100,              # Number of simulation iterations
  K = 1,                    # Number of source datasets (1 or more)
  n_source = 500,           # Number of source samples per source
  n_target_train = 100,     # Number of target training samples
  n_target_test = 100,      # Number of target test samples
  p1 = 15,                  # Number of source features
  p2 = 8,                   # Number of target features  
  regime = 'additive'       # Data generation regime
)
```

### K Parameter

- **K=1** (default): Single source scenario
  - Methods: HTL_Linear, HTL_Sieve, HmTL, Lasso_Target_Only
  - Generates 1 source dataset per iteration

- **K>1**: Multiple heterogeneous sources scenario
  - Methods: HTL_Linear_Multi, HTL_Sieve_Multi, HmTL_Multi, Lasso_Target_Only
  - Generates K independent source datasets per iteration
  - Aggregates information from all K sources

### Regime Options

- **'additive'**: Linear combination of first 5 features
  - Data: Y = X1 + X2 + X3 + X4 + X5 + noise

- **'multiplicative'**: Non-linear interactions
  - Data: Y = X1 * X2 + X3 + noise

- **'nonlinear'**: Complex non-linear patterns
  - Data: Y = sin(X1) + cos(X2) + X3² + noise

## Usage Examples

**Note:** Run all commands from the `Simulations/` directory.

### Option 1: Interactive R Session

```r
# Set working directory to Simulations folder first
setwd('./Simulations')

# Load and run simulation
source('run_simulation.R')

# Run with default parameters (50 iterations, additive regime)
results <- run_simulation(
  n_sim = 50,
  p1 = 15,
  p2 = 8,
  regime = 'additive',
  save_file = 'my_results.rds'
)

# Access results
head(results$summary)           # Individual prediction results
results$summary_stats           # Summary statistics by method
```

### Option 2: Command Line (Quick Run)

```bash
# Navigate to Simulations directory
cd Simulations

# Run simulation and save results
R --slave -e "
  source('simulation.R')
  set.seed(2026)
  results <- main(50, 500, 100, 100, 15, 8, 'additive')
  saveRDS(results, 'sim_results.rds')
"

# Summarize results
Rscript summarize_sim.R sim_results.rds

# Create plots
Rscript plot_sim_results.R sim_results.rds sim_plot.pdf
```

### Option 3: Complete Pipeline

```bash
# From Simulations/ directory
R -e "
  source('simulation.R')
  set.seed(2026)
  results <- main(100, 500, 100, 100, 15, 8, 'additive')
  saveRDS(results, 'sim_results.rds')
" && \
Rscript plot_sim_results.R sim_results.rds && \
Rscript summarize_sim.R sim_results.rds
```

## Output Interpretation

### Summary Statistics Table

The `summary_stats` data frame shows:
- **mean_rmse**: Average prediction error across test sets
- **sd_rmse**: Standard deviation of RMSE
- **median_rmse**: Median RMSE (robust to outliers)
- **mean_mae**: Mean absolute error
- **sd_mae**: Standard deviation of MAE
- **n**: Number of simulation iterations

### Performance Comparison

Results are compared relative to the **Lasso-Target-Only** baseline:
- **Positive improvement %**: Method performs better than baseline
- **Negative improvement %**: Method performs worse than baseline
- **Significance codes**: *** p<0.001, ** p<0.01, * p<0.05

## Expected Results

Typical findings from transfer learning comparisons:

1. **HTL-Sieve** often shows best performance when:
   - Feature relationships are non-linear
   - Source sample size is large
   - Number of proxy features (p1 - p2) is substantial

2. **HTL-Linear** typically outperforms:
   - Lasso-Target-Only baseline (when feature relationships are simple)
   - HmTL (when features are truly heterogeneous)

3. **HmTL** works well when:
   - Features are actually homogeneous (contrary to the heterogeneous setting)
   - Sample sizes are large enough to benefit from source data

4. **Lasso-Target-Only**:
   - Baseline performance
   - Often limited by small target sample size alone

## Algorithm Details

### HTL - Linear Mapping

1. Select common features between domains
2. Learn linear regression maps: f_j(X_common) → X_proxy_j for each proxy feature
   - Uses OLS: `glmnet(..., alpha=1, lambda=0)`
3. Predict proxy features on target data using learned maps
4. Train lasso on combined features (target + imputed proxy)

### HTL - Sieve Mapping

1. Same feature selection as linear approach
2. Learn non-linear mappings using sieve basis functions:
   - Cosine basis (or 'cosine' type in `sieve_preprocess()`)
   - Basis size selected via cross-validation
3. L1 regularization on basis coefficients
4. Predict and combine with target features same as linear approach

### Performance Metrics

- **RMSE**: Root Mean Squared Error - sensitive to large errors
  - RMSE = √(mean((Y_true - Y_pred)²))
  
- **MAE**: Mean Absolute Error - robust to outliers
  - MAE = mean(|Y_true - Y_pred|)

## References

- Zhang, T., & Simon, N. (2023). "Regression in tensor product spaces by the method of sieves." _Electronic Journal of Statistics_, 17(2), 3660.
- Original Sieve package: https://github.com/terrytianyuzhang/Sieve

## Notes

- Results may vary based on random seed
- Large simulations (n_sim=1000+) can be computationally intensive
- Parallel processing used in `cv_fitH_from_prxy()` for efficiency
- Requires: glmnet, data.table, tidyverse, Sieve, gridExtra (for plotting)
