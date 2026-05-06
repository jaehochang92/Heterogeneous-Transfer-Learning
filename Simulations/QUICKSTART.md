# Quick Start Guide - Transfer Learning Simulation

## What Was Created

This package supports transfer learning simulations with both **single source** (K=1) and **multiple sources** (K>1):

### Single Source (K=1) - 4 Methods
1. **HTL-Linear**: Heterogeneous TL with linear feature mapping
2. **HTL-Sieve**: Heterogeneous TL with non-linear sieve estimation  
3. **HmTL**: Homogeneous TL (same features across domains)
4. **Lasso-Only**: Baseline using target data only

### Multiple Sources (K>1) - 4 Methods
1. **HTL-Linear-Multi**: Aggregates linear mappings from K sources
2. **HTL-Sieve-Multi**: Aggregates sieve-based mappings from K sources
3. **HmTL-Multi**: Pools K homogeneous source datasets
4. **Lasso-Only**: Baseline using target data only

## File Summary

| File | Purpose |
|------|---------|
| `simulation.R` | Core: data generation, all methods, supports K parameter |
| `run_simulation.R` | Wrapper function with K support |
| `test_simulation.R` | Tests both K=1 and K>1 scenarios |
| `sv-cv.R` | Sieve basis estimation (dependency) |
| `sv-prim.R` | Non-linear basis functions (dependency) |
| `summarize_sim.R` | Summary statistics and t-tests |
| `plot_sim_results.R` | Visualization plots |
| `SIMULATION_README.md` | Detailed documentation |

## The Key Difference: HTL vs HmTL

```
Single Source (K=1):          Homogeneous (HmTL):
Source: X_s ∈ ℝ^(p1)         Source: X_s ∈ ℝ^p
Target: X_t ∈ ℝ^(p2)         Target: X_t ∈ ℝ^p
where p1 > p2                where same features

Solution: Learn mapping      Solution: Combine source
f: X_s → X̂_missing          & target data directly
```

Multiple Sources: Repeat mapping for each source and aggregate.

## Quick Start: Run a Simulation

### Setup (Required First Run)

```r
# Navigate to the Simulations directory
setwd('./Simulations')  # or navigate in terminal: cd Simulations
```

### Single-Source Simulation

```r
source('simulation.R')

# Run 50 simulations (single source, K=1)
results <- main(
  n_sim = 50,           
  K = 1,                # single source (default)
  n_source = 500,       
  n_target_train = 100, 
  n_target_test = 100,  
  p1 = 15,              # source features
  p2 = 8,               # target features
  regime = 'additive'   
)

print(results$summary_stats)
```

### Multi-Source Simulation

```r
source('simulation.R')

# Run 50 simulations with K=4 heterogeneous sources
results <- main(
  n_sim = 50,           
  K = 4,                # 4 heterogeneous source datasets
  n_source = 500,       
  n_target_train = 100, 
  n_target_test = 100,  
  p1 = 15,              # source features per source
  p2 = 8,               # target features
  regime = 'additive'   
)

print(results$summary_stats)
```

### Option 2: Save and Summarize

```bash
# Navigate to Simulations directory first
cd Simulations

# Single-source simulation (K=1)
R --slave -e "source('simulation.R'); set.seed(2026); 
              results <- main(n_sim=100, K=1, n_source=500, 
                            n_target_train=100, n_target_test=100, 
                            p1=15, p2=8, regime='additive')
              saveRDS(results, 'single_source.rds')"

# Multi-source simulation (K=4)
R --slave -e "source('simulation.R'); set.seed(2026); 
              results <- main(n_sim=100, K=4, n_source=500, 
                            n_target_train=100, n_target_test=100, 
                            p1=15, p2=8, regime='additive')
              saveRDS(results, 'multi_source_k4.rds')"

# Print summary
Rscript summarize_sim.R multi_source_k4.rds

# Create plots
Rscript plot_sim_results.R multi_source_k4.rds
```

### Option 3: Using Convenient Wrapper

```r
source('run_simulation.R')

# Single source
results <- run_simulation(n_sim=50, K=1, p1=15, p2=8, regime='additive')

# Multiple sources
results_multi <- run_simulation(n_sim=50, K=4, p1=15, p2=8, 
                                regime='additive',
                                save_file='multi_k4_results.rds')
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
