library(data.table)
library(ggplot2)

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
  dt[, SweepValue := fifelse(
    SweepParam == "K",
    as.character(K),
    fifelse(
      SweepParam == "np",
      as.character(np),
      fifelse(
        SweepParam == "nt",
        as.character(nt),
        fifelse(
          SweepParam == "p1",
          as.character(p1),
          fifelse(SweepParam == "p2", as.character(p2), as.character(dgp_fmap))
        )
      )
    )
  )]
  dt[, SweepKey := paste(SweepParam, SweepValue, sep = "::")]
  
  ordered_keys <- unlist(lapply(unique(dt$SweepParam), function(param) {
    vals <- unique(dt[SweepParam == param, SweepValue])
    num_vals <- suppressWarnings(as.numeric(vals))
    if (all(!is.na(num_vals))) {
      vals <- as.character(sort(unique(num_vals)))
    } else {
      vals <- sort(unique(vals))
    }
    paste(param, vals, sep = "::")
  }), use.names = FALSE)
  
  dt[, SweepKey := factor(SweepKey, levels = unique(ordered_keys))]
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
  
  long_dt[, .(Mean = mean(Value, na.rm = TRUE), N = .N), by = .(Metric, Method, Sp, FMap, SweepParam, SweepValue, SweepKey)]
}

plot_metric <- function(plot_dt, metric_name) {
  metric_dt <- plot_dt[Metric == metric_name]
  ggplot(metric_dt,
         aes(
           x = SweepKey,
           y = Mean,
           color = Method,
           shape = Method,
           group = Method
         )) +
    geom_line(linewidth = 0.7, alpha = 0.9) +
    geom_point(size = 2, alpha = 0.9) +
    facet_grid(
      rows = vars(SweepParam, FMap),
      cols = vars(Sp),
      scales = "free_y",
      # space = "free_y",
      labeller = label_parsed
    ) +
    scale_x_discrete(
      labels = function(x)
        sub("^[^:]+::", "", x)
    ) +
    labs(
      x = "Swept Value",
      y = sprintf("Mean %s", metric_name),
      color = "Method",
      title = sprintf("%s by Sweep Parameter, fmap, and Sparsity Case", metric_name)
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

run_plot_pipeline <- function(results_dir = "results",
                              out_dir = "results/plots") {
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)
  
  raw_dt <- load_simulation_results(results_dir)
  raw_dt <- add_sweep_columns(raw_dt)
  plot_dt <- build_plot_table(raw_dt)
  plot_dt <- plot_dt %>% 
    mutate(Sp = case_when(
      Sp == "sp"  ~ "sparse~delta",
      Sp == "nsp" ~ "dense~delta",
      TRUE        ~ Sp
    ))
  
  p_est <- plot_metric(plot_dt, "Est")
  p_pe <- plot_metric(plot_dt, "PE")
  
  ggsave(
    file.path(out_dir, "summary_est.png"),
    p_est,
    width = 14,
    height = 10,
    dpi = 300
  )
  ggsave(
    file.path(out_dir, "summary_pe.png"),
    p_pe,
    width = 14,
    height = 10,
    dpi = 300
  )
  
  list(est_plot = p_est,
       pe_plot = p_pe,
       summary = plot_dt)
}

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else "results"
out_dir <- if (length(args) >= 2) args[[2]] else "results/plots"

res <- run_plot_pipeline(results_dir = results_dir, out_dir = out_dir)
print(res$est_plot)
print(res$pe_plot)