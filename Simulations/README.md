# Simulations

This directory contains all simulation code for comparing heterogeneous transfer learning (HTL) methods.

## Quick Start

1. **Read the docs first**: Start with [QUICKSTART.md](QUICKSTART.md) for a 5-minute overview
2. **Detailed guide**: See [SIMULATION_README.md](SIMULATION_README.md) for complete documentation
3. **Run a test**: Try `Rscript test_simulation.R` to verify everything works

## Files

- **simulation.R** - Core simulation code (4 transfer learning methods)
- **run_simulation.R** - Convenient wrapper with preset configurations  
- **summarize_sim.R** - Statistical summary and comparison
- **plot_sim_results.R** - Visualization and plotting
- **test_simulation.R** - Quick validation test

## Prerequisites

All dependencies (sv-cv.R, sv-prim.R) are in the parent directory and automatically sourced.

Required R packages:
- glmnet
- tidyverse
- data.table
- Sieve (for sieve-based methods)

## Example Usage

```r
# From R console (in this directory)
source('simulation.R')
results <- main(n_sim=50, p1=15, p2=8, regime='additive')
print(results$summary_stats)
```

```bash
# From terminal (in this directory)
Rscript test_simulation.R  # Quick 5-iteration test
```

See QUICKSTART.md for more examples.
