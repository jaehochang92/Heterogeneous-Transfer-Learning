#!/usr/bin/env Rscript
# Quick test of simulation code - both single and multi-source

source('../sv-cv.R')
source('../sv-prim.R')
source('simulation.R')

cat("Testing Heterogeneous Transfer Learning Simulation\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Test 1: Single-source simulation (K=1)
cat("Test 1: Single-source simulation (K=1, 5 iterations)\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
set.seed(123)
results_single <- main(n_sim = 5, 
                       K = 1,
                       n_source = 200, 
                       n_target_train = 50, 
                       n_target_test = 50, 
                       p1 = 12, 
                       p2 = 6, 
                       regime = 'additive')

cat('\nSingle-Source Results:\n')
print(results_single$summary_stats)
cat("\n")

# Test 2: Multi-source simulation (K=3)
cat("Test 2: Multi-source simulation (K=3, 5 iterations)\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
set.seed(456)
results_multi <- main(n_sim = 5,
                      K = 3,
                      n_source = 200,
                      n_target_train = 50,
                      n_target_test = 50,
                      p1 = 12,
                      p2 = 6,
                      regime = 'additive')

cat('\nMulti-Source (K=3) Results:\n')
print(results_multi$summary_stats)
cat("\n")

cat("Tests completed successfully!\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
