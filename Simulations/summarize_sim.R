# Summary Script for Transfer Learning Simulation
# Arguments: takes RDS file path as input

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript summarize_sim.R <results_file.rds>")
}

results_file <- args[1]

if (!file.exists(results_file)) {
  stop(paste("File not found:", results_file))
}

library(data.table)

results <- readRDS(results_file)

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SIMULATION RESULTS SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Print summary statistics
cat("Summary Statistics by Method:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
print(results$summary_stats)
cat("\n")

# Calculate improvements relative to baseline (Lasso Target Only)
baseline_rmse <- results$summary_stats[method == "Lasso_Target_Only", mean_rmse]
baseline_mae <- results$summary_stats[method == "Lasso_Target_Only", mean_mae]

results$summary_stats[, 
                      rmse_improvement_pct := ((baseline_rmse - mean_rmse) / baseline_rmse) * 100]
results$summary_stats[, 
                      mae_improvement_pct := ((baseline_mae - mean_mae) / baseline_mae) * 100]

cat("Performance Improvement (%) over Lasso Target Only Baseline:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
improvement <- results$summary_stats[, .(method, rmse_improvement_pct, mae_improvement_pct)]
colnames(improvement) <- c("Method", "RMSE Improvement (%)", "MAE Improvement (%)")
print(improvement)
cat("\n")

# Detailed method comparison
cat("Detailed Method Comparison:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
for (method_name in unique(results$summary_stats$method)) {
  stats <- results$summary_stats[method == method_name]
  cat(sprintf("%s:\n", method_name))
  cat(sprintf("  RMSE: %.4f ± %.4f (Median: %.4f)\n",
              stats$mean_rmse, stats$sd_rmse, stats$median_rmse))
  cat(sprintf("  MAE:  %.4f ± %.4f (Median: %.4f)\n",
              stats$mean_mae, stats$sd_mae, stats$median_mae))
  cat(sprintf("  N samples: %d\n\n", stats$n))
}

# Statistical significance tests (pairwise t-tests)
cat("Pairwise T-Tests (vs Lasso Target Only):\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

baseline_rmses <- results$summary[method == "Lasso_Target_Only", rmse]

for (method_name in unique(results$summary_stats$method)) {
  if (method_name != "Lasso_Target_Only") {
    method_rmses <- results$summary[method == method_name, rmse]
    t_test <- t.test(method_rmses, baseline_rmses, paired = TRUE)
    sig <- ifelse(t_test$p.value < 0.001, "***",
                  ifelse(t_test$p.value < 0.01, "**",
                         ifelse(t_test$p.value < 0.05, "*", "ns")))
    cat(sprintf("  %s vs Lasso_Target_Only: t = %7.4f, p-value = %.4e %s\n",
                method_name, t_test$statistic, t_test$p.value, sig))
  }
}
cat("\nSignificance codes: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Simulation completed successfully.\n")
cat(sprintf("Total simulations: %d\n", nrow(results$summary) / length(unique(results$summary$method))))
cat(paste(rep("=", 80), collapse = ""), "\n\n")
