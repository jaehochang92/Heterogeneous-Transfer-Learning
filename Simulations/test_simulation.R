#!/usr/bin/env Rscript
# Quick test of simulation code

source('../sv-cv.R')
source('../sv-prim.R')
source('simulation.R')

set.seed(123)
cat('Running quick test with 5 simulations...\n')
results <- main(n_sim = 5, 
                n_source = 200, 
                n_target_train = 50, 
                n_target_test = 50, 
                p1 = 12, 
                p2 = 6, 
                regime = 'additive')

cat('\nTest Summary:\n')
print(results$summary_stats)

cat('\nTest completed successfully!\n')
cat('Result object contains:\n')
cat('  $results: List of individual simulation results\n')
cat('  $summary: Data table with RMSE and MAE for each iteration\n')
cat('  $summary_stats: Summary statistics by method\n')
