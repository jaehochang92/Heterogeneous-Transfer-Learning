# Quick Start Guide - Transfer Learning Simulation

## What Was Created

This package implements a comprehensive simulation comparing **4 methods** for transfer learning:

1. **HTL-Linear**: Heterogeneous TL with linear feature mapping
2. **HTL-Sieve**: Heterogeneous TL with non-linear sieve estimation  
3. **HmTL**: Homogeneous TL (same features across domains)
4. **Lasso-Only**: Baseline using target data only

## File Summary

| File | Purpose |
|------|---------|
| `simulation.R` | Core simulation: data generation, all 4 methods, main loop |
| `sv-cv.R` | Sieve basis estimation with cross-validation (dependency) |
| `sv-prim.R` | Non-linear basis functions (dependency) |
| `run_simulation.R` | Convenient wrapper to launch simulations |
| `summarize_sim.R` | Print summary statistics and t-tests |
| `plot_sim_results.R` | Create visualization plots |
| `test_simulation.R` | Quick test to verify everything works |
| `SIMULATION_README.md` | Detailed documentation |

## The Key Difference: HTL vs HmTL

```
Heterogeneous (HTL):          Homogeneous (HmTL):
Source: X_s ∈ ℝ^(p1)         Source: X_s ∈ ℝ^p
Target: X_t ∈ ℝ^(p2)         Target: X_t ∈ ℝ^p
where p1 > p2                where same features

Solution: Learn mapping      Solution: Combine source
f: X_s → X̂_missing          & target data directly
```

## Quick Start: Run a Simulation

### Setup (Required First Run)

```r
# Navigate to the Simulations directory
setwd('./Simulations')  # or navigate in terminal: cd Simulations
```

### Option 1: Simple Interactive Run

```r
source('simulation.R')

# Run 50 simulations (5-10 minutes on modern hardware)
results <- main(
  n_sim = 50,           
  n_source = 500,       
  n_target_train = 100, 
  n_target_test = 100,  
  p1 = 15,              # source features
  p2 = 8,               # target features
  regime = 'additive'   
)

# View results
results$summary_stats  # Summary table
```

### Option 2: Save and Summarize (from Simulations/ directory)

```bash
# Navigate to Simulations directory first
cd Simulations

# Run simulation and save results
R --slave -e "
  source('simulation.R')
  set.seed(2026)
  results <- main(100, 500, 100, 100, 15, 8, 'additive')
  saveRDS(results, 'my_results.rds')
"

# Print summary
Rscript summarize_sim.R my_results.rds

# Create plots
Rscript plot_sim_results.R my_results.rds my_plots.pdf
```

### Option 3: Using Convenient Wrapper

```r
source('run_simulation.R')

results <- run_simulation(
  n_sim = 50,
  regime = 'nonlinear',
  save_file = 'htl_vs_baseline.rds'
)
```

## Understanding the Output

### Summary Statistics Table:
```
           method mean_rmse  sd_rmse median_rmse mean_mae  sd_mae median_mae  n
1    HTL_Linear      0.523   0.115       0.502    0.398   0.089       0.382 50
2     HTL_Sieve      0.498   0.108       0.489    0.379   0.085       0.367 50
3          HmTL      0.487   0.102       0.476    0.372   0.080       0.361 50
4 Lasso_Target_Only 0.615   0.142       0.601    0.468   0.108       0.453 50
```

### Key Metrics:
- **RMSE**: Root Mean Squared Error (prediction accuracy)
- **MAE**: Mean Absolute Error (robust to outliers)
- Lower values are better

### Interpreting Results:
- **HTL-Sieve RMSE 0.498** vs **Lasso-Only 0.615** = **19% improvement**
- Statistical significance shows if improvement is just noise

## Customization Options

### Different Regime Types:
```r
# Linear relationship (easiest to learn)
main(n_sim=50, regime='additive')     

# Multiplicative interactions  
main(n_sim=50, regime='multiplicative')

# Non-linear transformations (hardest to learn)
main(n_sim=50, regime='nonlinear')
```

### Different Data Sizes:
```r
# Small sample test (fast, ~2-3 min)
main(n_sim=10, n_source=200, n_target_train=50, p1=10, p2=5)

# Standard simulation (medium, ~5-10 min)
main(n_sim=50, n_source=500, n_target_train=100, p1=15, p2=8)  

# Large simulation (slow, ~30-60 min)
main(n_sim=100, n_source=1000, n_target_train=150, p1=20, p2=10)
```

## Expected Results Summary

### Typical Performance Order (Best → Worst):
1. **HTL-Sieve** - Best when relationships are non-linear
2. **HmTL** - Good when features are actually shared
3. **HTL-Linear** - Good for simple relationships
4. **Lasso-Only** - Baseline (limited by small target sample)

### When Each Method Excels:

| Method | Best When | Advantage |
|--------|-----------|-----------|
| HTL-Sieve | Complex non-linear feature relationships | Captures intricate patterns |
| HmTL | Features truly align across domains | Simple, no feature mapping needed |
| HTL-Linear | Linear feature relationships | Computationally efficient |
| Lasso-Only | When transfer learning fails | Fast, no complex mappings |

## Troubleshooting

### "Error: could not find function X"
→ Make sure you ran `source('simulation.R')` after `source('sv-cv.R')`

### "Sieve package not found"
→ Install: `install.packages('Sieve')`

### Simulation runs very slow
→ Reduce n_sim (try 10), n_source (try 200), or both

### Poor results for all methods
→ May mean simulated data is too simple; try `regime='nonlinear'`

## Advanced: Parallel Processing

The feature mapping (`cv_fitH_from_prxy`) uses parallel processing:
- Automatically detects CPU cores
- Uses `pbmcapply` if available (shows progress bar)
- Falls back to sequential if needed

Install for progress tracking:
```r
install.packages('pbmcapply')
```

## References

- Based on: Zhang, T., & Simon, N. (2023). Regression in tensor product spaces by the method of sieves.
- Sieve package: https://github.com/terrytianyuzhang/Sieve

## Contact

For questions about model implementation or results interpretation, refer to:
- `SIMULATION_README.md` - Detailed documentation
- `simulation.R` - Heavily commented source code
- Method references in script headers
