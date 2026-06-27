library(dplyr)
library(data.table)
library(ggplot2)

log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

resolve_sweep_value <- function(dt) {
  result <- as.character(dt$dgp_fmap)
  result[dt$SweepParam == "h"] <- as.character(dt$h[dt$SweepParam == "h"])
  result[dt$SweepParam == "K"] <- as.character(dt$K[dt$SweepParam == "K"])
  result[dt$SweepParam == "np"] <- as.character(dt$np[dt$SweepParam == "np"])
  result[dt$SweepParam == "nt"] <- as.character(dt$nt[dt$SweepParam == "nt"])
  result[dt$SweepParam == "p1"] <- as.character(dt$p1[dt$SweepParam == "p1"])
  result[dt$SweepParam == "p2"] <- as.character(dt$p2[dt$SweepParam == "p2"])
  result
}

build_sweep_key_levels <- function(dt) {
  ordered_keys <- character(0)
  for (param in unique(dt$SweepParam)) {
    values <- unique(dt[dt$SweepParam == param, ]$SweepValue)
    numeric_values <- suppressWarnings(as.numeric(values))
    if (all(!is.na(numeric_values))) {
      values <- as.character(sort(unique(numeric_values)))
    } else {
      values <- sort(unique(values))
    }
    ordered_keys <- c(ordered_keys, paste(param, values, sep = "::"))
  }
  unique(ordered_keys)
}

summarize_plot_values <- function(long_dt) {
  long_dt <- copy(long_dt)
  long_dt[Metric %in% c("Est", "PE") & Value > 0, Value := log(Value)]
  # long_dt[Metric %in% c("Est", "PE") & Value <= 0, Value := NA]
  
  summary_dt <- long_dt[, .(
    Mean = mean(Value, na.rm = TRUE),
    LQ   = quantile(Value, 0.25, na.rm = TRUE),
    UQ   = quantile(Value, 0.75, na.rm = TRUE),
    N    = sum(!is.na(Value))
  ), by = .(Metric, Method, Sp, FMap, SweepParam, SweepValue, SweepKey)]
  
  return(summary_dt)
}

load_simulation_results <- function(results_dir = "results_sand") {
  batch_dir <- file.path(results_dir, "batch_bundles")
  if (!dir.exists(batch_dir)) {
    stop(sprintf("ERROR: dir '%s' does not exist.", batch_dir))
  }
  
  files <- list.files(batch_dir, pattern = "^sand_res_(bundle|batch)_.*[.]rds$", full.names = TRUE)
  files <- files[!grepl("\\.partial\\.rds$", files)]
  
  if (length(files) == 0L) {
    stop(sprintf("ERROR: dir '%s' has no *.rds files.", batch_dir))
  }
  
  log_info(sprintf("Aggregating %d batch data files...", length(files)))
  
  rbindlist(lapply(files, readRDS), fill = TRUE, use.names = TRUE)
}

add_sweep_columns <- function(dt) {
  dt <- copy(dt)
  dt[, SweepValue := resolve_sweep_value(dt)]
  dt[, SweepKey := paste(SweepParam, SweepValue, sep = "::")]
  dt[, SweepKey := factor(SweepKey, levels = build_sweep_key_levels(dt))]
  dt
}

build_plot_table <- function(dt) {
  metric_cols <- names(dt)[grepl("^(Est|PE)_|Precision|Recall|F1|ExactMatch|FalsePos|FalseNeg|N_selected", names(dt))]
  if (length(metric_cols) == 0) {
    stop("No Est_/PE_ columns found in results.")
  }
  
  if (!("dgp_fmap" %in% names(dt))) {
    stop("Column 'dgp_fmap' is required to split plots by fmap.")
  }
  long_dt <- melt(
    dt,
    id.vars = c("Sp", "dgp_fmap", "SweepParam", "SweepValue", "SweepKey"),
    measure.vars = metric_cols,
    variable.name = "MetricMethod",
    value.name = "Value"
  )
  
  long_dt[, Metric := sub("_.*", "", MetricMethod)]
  long_dt[, Method := sub("^[^_]+_", "", MetricMethod)]
  long_dt[is.na(Method), Method := "SAND"]
  
  long_dt[, FMap := as.character(dgp_fmap)]
  summarize_plot_values(long_dt)
}

plot_metric <- function(plot_dt, metric_name, sweep_param = NULL, fmap_val = NULL) {
  metric_dt <- plot_dt[plot_dt$Metric == metric_name]
  if (!is.null(sweep_param)) {
    metric_dt <- metric_dt[metric_dt$SweepParam == sweep_param]
  }
  if (is.factor(metric_dt$SweepKey)) {
    metric_dt[, SweepKey := droplevels(SweepKey)]
  }
  if (!is.null(fmap_val)) {
    metric_dt <- metric_dt[FMap == fmap_val]
    if (fmap_val == "linear") {
      metric_dt <- metric_dt[Method != "HTL.NL"]
    } else if (fmap_val == "nlinear") {
      metric_dt <- metric_dt[Method != "HTL.LN"]
    }
  }
  ggplot(metric_dt,
         aes(
           x = SweepKey,
           y = Mean,
           color = Method,
           shape = Method,
           group = interaction(Method, Sp)
         )) +
         geom_errorbar(
      aes(ymin = LQ, ymax = UQ),
      width = 0.1,
      alpha = 0.35,
      linewidth = 0.5
    ) +
    geom_line(linewidth = 1, alpha = 0.5) +
    geom_point(size = 4, alpha = 0.5) +
    facet_grid(
      cols = vars(Sp),
      scales = "free_y",
      labeller = label_parsed
    ) +
    scale_x_discrete(
      labels = function(x)
        sub("^[^:]+::", "", x)
    ) +
    labs(
      x = if (is.null(sweep_param)) "Swept Value" else sweep_param,
      y = ifelse(metric_name == 'Est', 'LogRMS(EE)', 'LogRMS(PE)'),
      color = "Method",
      shape = "Method",
      title = if (is.null(sweep_param)) {
        sprintf("%s", metric_name)
      } else {
        sprintf("%s for growing %s", ifelse(
          metric_name == 'Est', 
          'Estimation Errors (EE)', 'Prediction Errors (PE)'), 
          sweep_param)
      }
    ) +
    theme_minimal(base_size = 30) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = 'bottom'
    )
}

plot_selection_metrics <- function(plot_dt, sweep_param = NULL, fmap_val = NULL) {
  sel_dt <- plot_dt[Metric %in% c("Precision", "Recall", "F1", "ExactMatch", "FalsePos")]
  
  if (!is.null(sweep_param)) sel_dt <- sel_dt[SweepParam == sweep_param]
  if (is.factor(sel_dt$SweepKey)) {
    sel_dt[, SweepKey := droplevels(SweepKey)]
  }
  if (!is.null(fmap_val)) sel_dt <- sel_dt[FMap == fmap_val]
  
  ggplot(sel_dt, aes(x = SweepKey, y = Mean, color = Metric, group = Metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    facet_grid(Metric ~ Sp, scales = "free_y", labeller = label_parsed) +
    labs(
      x = ifelse(is.null(sweep_param), "Swept Value", sweep_param),
      y = "Performance Score",
      title = "Source Selection Accuracy"
    ) +
    theme_minimal(base_size = 30) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
}

run_plot_pipeline <- function(results_dir = "results_sand", 
                              out_dir = "results_sand/plots",
                              exclude = NULL) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  raw_dt <- load_simulation_results(results_dir)
  raw_dt <- add_sweep_columns(raw_dt)
  if (!is.null(exclude)) {
    plot_dt <- raw_dt %>% dplyr::select(!contains(exclude)) %>% build_plot_table
  } else plot_dt <- raw_dt %>% build_plot_table
  plot_dt[Sp == "sp", Sp := "sparse~delta"]
  plot_dt[Sp == "nsp", Sp := "dense~delta"]
  
  sweep_params <- unique(as.character(plot_dt$SweepParam))
  fmap_cases <- unique(as.character(plot_dt$FMap))
  est_plots <- vector("list", length(sweep_params))
  pe_plots <- vector("list", length(sweep_params))
  names(est_plots) <- sweep_params
  names(pe_plots) <- sweep_params
  
  for (sp_name in sweep_params) {
    for (fmap_val in fmap_cases) {
      plot_key <- paste(sp_name, fmap_val, sep = "_")
      
      p_est <- plot_metric(plot_dt, "Est", sweep_param = sp_name, 
                           fmap_val = fmap_val)
      p_pe <- plot_metric(plot_dt, "PE", sweep_param = sp_name, 
                          fmap_val = fmap_val)
      
      est_plots[[plot_key]] <- p_est
      pe_plots[[plot_key]] <- p_pe
      
      ggsave(
        file.path(out_dir, sprintf("ee_%s_%s.pdf", sp_name, fmap_val)),
        p_est,
        width = 12,
        height = 8,
        dpi = 300
      )
      ggsave(
        file.path(out_dir, sprintf("pe_%s_%s.pdf", sp_name, fmap_val)),
        p_pe,
        width = 12,
        height = 8,
        dpi = 300
      )
      p_sel <- plot_selection_metrics(plot_dt, sweep_param = sp_name, fmap_val = fmap_val)
      ggsave(
        file.path(out_dir, sprintf("sel_%s_%s.pdf", sp_name, fmap_val)),
        p_sel, width = 12, height = 10, dpi = 300
      )
    }
  }
  
  list(
    est_plots = est_plots,
    pe_plots = pe_plots,
    summary = plot_dt
  )
}

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else "results_sand"
out_dir <- if (length(args) >= 2) args[[2]] else file.path(results_dir, "plots")

res <- run_plot_pipeline(results_dir = results_dir, out_dir = out_dir, 
                         exclude = c('OracleHmTL'))