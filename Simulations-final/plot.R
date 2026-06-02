library(data.table)
library(ggplot2)

resolve_sweep_value <- function(dt) {
  result <- as.character(dt$dgp_fmap)
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
  summary_dt <- long_dt[, Value := log(Value + 1)]
  summary_dt <- long_dt[, .(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    N = .N
  ), by = .(Metric, Method, Sp, FMap, SweepParam, SweepValue, SweepKey)]
  summary_dt[, SE := SD / sqrt(N)]
  summary_dt
}

load_simulation_results <- function(results_dir = "results") {
  files <- list.files(results_dir, pattern = "^sim_res_task_[0-9]{4}[.]rds$", full.names = TRUE)
  if (length(files) == 0) {
    stop(sprintf("No result files found in '%s'.", results_dir))
  }
  rbindlist(lapply(files, readRDS),
            fill = TRUE,
            use.names = TRUE)
}

add_sweep_columns <- function(dt) {
  dt <- copy(dt)
  dt[, SweepValue := resolve_sweep_value(dt)]
  dt[, SweepKey := paste(SweepParam, SweepValue, sep = "::")]
  dt[, SweepKey := factor(SweepKey, levels = build_sweep_key_levels(dt))]
  dt
}

build_plot_table <- function(dt) {
  metric_cols <- grep("^(Est|PE)_", names(dt), value = TRUE)
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
    value.name = "Value",
    variable.factor = FALSE
  )
  long_dt[, Metric := sub("_.*", "", MetricMethod)]
  long_dt[, Method := sub("^[^_]+_", "", MetricMethod)]
  long_dt[, FMap := as.character(dgp_fmap)]
  
  summarize_plot_values(long_dt)
}

plot_metric <- function(plot_dt, metric_name, sweep_param = NULL) {
  metric_dt <- plot_dt[plot_dt$Metric == metric_name]
  if (!is.null(sweep_param)) {
    metric_dt <- metric_dt[metric_dt$SweepParam == sweep_param]
  }
  ggplot(metric_dt,
         aes(
           x = SweepKey,
           y = Mean,
           color = Method,
           shape = Method,
           group = Method
         )) +
    geom_errorbar(
      aes(ymin = Mean - 2 * SE, ymax = Mean + 2 * SE),
      width = 0.15,
      alpha = 0.45,
      linewidth = 0.5
    ) +
    geom_line(linewidth = 1, alpha = 0.9) +
    geom_point(size = 4, alpha = 0.9) +
    facet_grid(
      rows = vars(FMap),
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
      y = sprintf("Mean log(%s)", metric_name),
      color = "Method",
      title = if (is.null(sweep_param)) {
        sprintf("%s", metric_name)
      } else {
        sprintf("%s for sweep %s", metric_name, sweep_param)
      }
    ) +
    theme_bw(base_size = 18) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

run_plot_pipeline <- function(results_dir = "results",
                              out_dir = "results/plots") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  raw_dt <- load_simulation_results(results_dir)
  raw_dt <- add_sweep_columns(raw_dt)
  plot_dt <- build_plot_table(raw_dt)
  plot_dt[Sp == "sp", Sp := "sparse~delta"]
  plot_dt[Sp == "nsp", Sp := "dense~delta"]

  sweep_params <- unique(as.character(plot_dt$SweepParam))
  est_plots <- vector("list", length(sweep_params))
  pe_plots <- vector("list", length(sweep_params))
  names(est_plots) <- sweep_params
  names(pe_plots) <- sweep_params

  for (sp_name in sweep_params) {
    p_est <- plot_metric(plot_dt, "Est", sweep_param = sp_name)
    p_pe <- plot_metric(plot_dt, "PE", sweep_param = sp_name)

    est_plots[[sp_name]] <- p_est
    pe_plots[[sp_name]] <- p_pe

    ggsave(
      file.path(out_dir, sprintf("summary_est_%s.png", sp_name)),
      p_est,
      width = 12,
      height = 8,
      dpi = 300
    )
    ggsave(
      file.path(out_dir, sprintf("summary_pe_%s.png", sp_name)),
      p_pe,
      width = 12,
      height = 8,
      dpi = 300
    )
  }

  list(
    est_plots = est_plots,
    pe_plots = pe_plots,
    summary = plot_dt
  )
}

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else "results"
out_dir <- if (length(args) >= 2) args[[2]] else "results/plots"

res <- run_plot_pipeline(results_dir = results_dir, out_dir = out_dir)