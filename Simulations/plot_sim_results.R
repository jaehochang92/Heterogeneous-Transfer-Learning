# Visualization and Summary of Simulation Results
# Compares performance metrics across all methods

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)

# Function to load and summarize results
summarize_results <- function(results_rds_file) {
  results <- readRDS(results_rds_file)
  
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Simulation Results Summary\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  # Print summary statistics
  cat("Summary Statistics by Method:\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  print(results$summary_stats)
  cat("\n")
  
  # Calculate improvements relative to baseline (Lasso Target Only)
  baseline_rmse <- results$summary_stats[method == "Lasso_Target_Only", mean_rmse]
  results$summary_stats[, 
                        rmse_improvement := ((baseline_rmse - mean_rmse) / baseline_rmse) * 100]
  
  cat("RMSE Improvement (%) over Lasso Target Only baseline:\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  improvement <- results$summary_stats[, .(method, rmse_improvement)]
  print(improvement)
  cat("\n")
  
  # Create detailed summary
  cat("Detailed Method Comparison:\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  for (method_name in unique(results$summary_stats$method)) {
    stats <- results$summary_stats[method == method_name]
    cat(sprintf(
      "%s:\n",
      method_name
    ))
    cat(sprintf(
      "  RMSE: %.4f (SD: %.4f) [Median: %.4f]\n",
      stats$mean_rmse, stats$sd_rmse, stats$median_rmse
    ))
    cat(sprintf(
      "  MAE:  %.4f (SD: %.4f) [Median: %.4f]\n",
      stats$mean_mae, stats$sd_mae, stats$median_mae
    ))
    cat("\n")
  }
  
  # Statistical significance test (t-tests)
  cat("Pairwise T-Tests (vs Lasso Target Only):\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  baseline_rmses <- results$summary[method == "Lasso_Target_Only", rmse]
  
  for (method_name in unique(results$summary_stats$method)) {
    if (method_name != "Lasso_Target_Only") {
      method_rmses <- results$summary[method == method_name, rmse]
      t_test <- t.test(method_rmses, baseline_rmses, paired = TRUE)
      cat(sprintf(
        "%s vs Lasso_Target_Only: t = %.4f, p-value = %.4e\n",
        method_name, t_test$statistic, t_test$p.value
      ))
    }
  }
  cat("\n")
  
  return(results)
}

# Plotting function
plot_results <- function(results, output_file = "sim_results_plot.pdf") {
  # Prepare data for plotting
  plot_dt <- results$summary
  
  # RMSE boxplot
  p1 <- ggplot(plot_dt, aes(x = method, y = rmse, fill = method)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    labs(
      title = "RMSE Distribution by Method",
      x = "Method",
      y = "RMSE (Test Set)",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # MAE boxplot
  p2 <- ggplot(plot_dt, aes(x = method, y = mae, fill = method)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    labs(
      title = "MAE Distribution by Method",
      x = "Method",
      y = "MAE (Test Set)",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Mean performance comparison
  summary_stats <- results$summary_stats
  p3 <- ggplot(summary_stats, aes(x = method, y = mean_rmse, fill = method)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_errorbar(
      aes(ymin = mean_rmse - sd_rmse, ymax = mean_rmse + sd_rmse),
      width = 0.2
    ) +
    labs(
      title = "Mean RMSE with Standard Deviation",
      x = "Method",
      y = "Mean RMSE",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Violin plot comparing all methods
  p4 <- ggplot(plot_dt, aes(x = method, y = rmse, fill = method)) +
    geom_violin(alpha = 0.7) +
    labs(
      title = "RMSE Distribution (Violin Plot)",
      x = "Method",
      y = "RMSE",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Combine plots
  combined_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  # Save
  if (!is.null(output_file)) {
    ggsave(output_file, combined_plot, width = 12, height = 10)
    cat(sprintf("Plot saved to %s\n", output_file))
  }
  
  return(list(p_rmse_box = p1, p_mae_box = p2, p_rmse_bar = p3, p_rmse_violin = p4))
}

# Main execution
if (interactive()) {
  # When running interactively, provide instructions
  cat("Load results using: results <- summarize_results('sim_results.rds')\n")
  cat("Create plots using: plots <- plot_results(results, 'results_plot.pdf')\n")
} else {
  # When called from command line
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    cat("Usage: Rscript plot_sim_results.R <results_file.rds> [output_file.pdf]\n")
    stop("Missing required argument: results file path")
  }
  
  results_file <- args[1]
  output_file <- if (length(args) > 1) args[2] else "sim_results_plot.pdf"
  
  if (!file.exists(results_file)) {
    stop(paste("File not found:", results_file))
  }
  
  results <- summarize_results(results_file)
  plots <- plot_results(results, output_file)
}
