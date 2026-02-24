# Main Runner Script for Transfer Learning Simulation
# Sources all required files and runs the complete analysis pipeline

cat("Loading libraries and source files...\n")

library(glmnet)
library(tidyverse)
library(data.table)

# Source required files
source('../sv-cv.R')
source('../sv-prim.R')
source('simulation.R')

# Main execution function
run_simulation <- function(n_sim = 50,
                          n_source = 500,
                          n_target_train = 100,
                          n_target_test = 100,
                          p1 = 15,
                          p2 = 8,
                          regime = 'additive',
                          save_file = 'sim_results.rds',
                          verbose = TRUE) {
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("TRANSFER LEARNING SIMULATION\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
  }
  
  # Run simulation
  results <- main(
    n_sim = n_sim,
    n_source = n_source,
    n_target_train = n_target_train,
    n_target_test = n_target_test,
    p1 = p1,
    p2 = p2,
    regime = regime
  )
  
  # Save results
  saveRDS(results, save_file)
  if (verbose) {
    cat(sprintf("Results saved to: %s\n\n", save_file))
  }
  
  return(results)
}

# Default execution (main function is already sourced from simulation.R)
# Uncomment one of the following to run with different configurations:

# Configuration 1: Small simulation for quick testing
# results <- run_simulation(
#   n_sim = 10,
#   n_source = 300,
#   n_target_train = 50,
#   n_target_test = 50,
#   p1 = 10,
#   p2 = 5,
#   regime = 'additive'
# )

# Configuration 2: Medium simulation (default)
# results <- run_simulation(
#   n_sim = 50,
#   n_source = 500,
#   n_target_train = 100,
#   n_target_test = 100,
#   p1 = 15,
#   p2 = 8,
#   regime = 'additive'
# )

# Configuration 3: Large simulation
# results <- run_simulation(
#   n_sim = 100,
#   n_source = 500,
#   n_target_train = 150,
#   n_target_test = 150,
#   p1 = 20,
#   p2 = 10,
#   regime = 'nonlinear'
# )

cat("Simulation runner loaded. Available functions:\n")
cat("  main() - Run simulation from simulation.R\n")
cat("  run_simulation() - Run with parameter control and save results\n\n")
cat("Example usage:\n")
cat("  results <- run_simulation(n_sim=50, p1=15, p2=8, regime='additive')\n")
