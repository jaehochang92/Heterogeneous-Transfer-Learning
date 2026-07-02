# Aggregate pilot CV results and compare SAND vs TLGLM source selection.
#
# Expects cache/resdf_<bundle_id>.rds from 2.main.R with columns:
#   sand_sel / tlglm_sel  (list columns of integer source indices), or
#   SAND_hA / TLGLM_hA    (legacy list columns)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

pilot_dir <- Sys.getenv(
  "SLURM_SUBMIT_DIR",
  unset = normalizePath(".", winslash = "/")
)
setwd(pilot_dir)

load_pilot_results <- function(cache_dir = "cache", pattern = "^resdf_.*\\.rds$") {
  files <- list.files(cache_dir, pattern = pattern, full.names = TRUE)
  files <- files[!grepl("resdf_all\\.rds$", files)]

  if (!length(files)) {
    stop(
      "No result files found in ", cache_dir,
      " matching ", pattern, ". Run sbmt.sh / 2.main.R first.",
      call. = FALSE
    )
  }

  bind_rows(lapply(files, function(f) {
    x <- readRDS(f)
    x$source_file <- basename(f)
    x
  }))
}

as_sel_vec <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(integer())
  }
  if (is.list(x) && length(x) == 1L) {
    x <- x[[1L]]
  }
  if (identical(x, "all")) {
    return(NA_integer_)
  }
  sort(unique(as.integer(x)))
}

pick_sel_column <- function(df, primary, legacy) {
  if (primary %in% names(df)) {
    return(primary)
  }
  if (legacy %in% names(df)) {
    return(legacy)
  }
  NA_character_
}

enrich_selection <- function(wide_df, K = NULL) {
  sand_col <- pick_sel_column(wide_df, "sand_sel", "SAND_hA")
  tlglm_col <- pick_sel_column(wide_df, "tlglm_sel", "TLGLM_hA")

  if (is.na(sand_col) || is.na(tlglm_col)) {
    stop(
      "Selection columns not found. Expected sand_sel/tlglm_sel ",
      "or SAND_hA/TLGLM_hA.",
      call. = FALSE
    )
  }

  out <- wide_df %>%
    mutate(
      sand_ids = lapply(.data[[sand_col]], as_sel_vec),
      tlglm_ids = lapply(.data[[tlglm_col]], as_sel_vec),
      sand_n = if ("sand_n" %in% names(.)) sand_n else vapply(sand_ids, length, integer(1)),
      tlglm_n = if ("tlglm_n" %in% names(.)) tlglm_n else vapply(tlglm_ids, length, integer(1)),
      jaccard = mapply(function(a, b) {
        if (length(a) == 0L && length(b) == 0L) {
          return(1)
        }
        if (length(a) == 0L || length(b) == 0L) {
          return(0)
        }
        length(intersect(a, b)) / length(union(a, b))
      }, sand_ids, tlglm_ids, SIMPLIFY = TRUE, USE.NAMES = FALSE),
      exact_match = mapply(function(a, b) {
        as.integer(identical(a, b))
      }, sand_ids, tlglm_ids, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    )

  if (!is.null(K)) {
    out <- out %>%
      mutate(
        sand_n = ifelse(vapply(sand_ids, function(x) length(x) == 1L && is.na(x[1]), logical(1)), K, sand_n),
        tlglm_n = ifelse(vapply(tlglm_ids, function(x) length(x) == 1L && is.na(x[1]), logical(1)), K, tlglm_n)
      )
  }

  out
}

infer_proxy_names <- function(proxy_path = "proxy.rds") {
  if (!file.exists(proxy_path)) {
    return(NULL)
  }
  nm <- names(readRDS(proxy_path))
  if (length(nm)) {
    return(nm)
  }
  NULL
}

selection_count_long <- function(sel_df) {
  sel_df %>%
    transmute(
      replicate,
      fold,
      cv_id = paste(replicate, fold, sep = "_"),
      SAND = sand_n,
      TLGLM = tlglm_n
    ) %>%
    pivot_longer(
      cols = c(SAND, TLGLM),
      names_to = "method",
      values_to = "n_selected"
    ) %>%
    mutate(
      method = factor(method, levels = c("SAND", "TLGLM"))
    )
}

selection_frequency <- function(sel_df, proxy_names = NULL) {
  K <- max(
    vapply(sel_df$sand_ids, function(x) if (length(x)) max(x) else 0L, integer(1)),
    vapply(sel_df$tlglm_ids, function(x) if (length(x)) max(x) else 0L, integer(1)),
    na.rm = TRUE
  )

  if (is.null(proxy_names) || length(proxy_names) < K) {
    proxy_names <- paste0("Source ", seq_len(K))
  }

  sand_mat <- matrix(0L, nrow = nrow(sel_df), ncol = K)
  tlglm_mat <- matrix(0L, nrow = nrow(sel_df), ncol = K)

  for (i in seq_len(nrow(sel_df))) {
    if (length(sel_df$sand_ids[[i]])) {
      sand_mat[i, sel_df$sand_ids[[i]]] <- 1L
    }
    if (length(sel_df$tlglm_ids[[i]])) {
      tlglm_mat[i, sel_df$tlglm_ids[[i]]] <- 1L
    }
  }

  tibble(
    source_id = seq_len(K),
    source = proxy_names[seq_len(K)],
    SAND = colMeans(sand_mat),
    TLGLM = colMeans(tlglm_mat)
  ) %>%
    pivot_longer(
      cols = c(SAND, TLGLM),
      names_to = "method",
      values_to = "selection_rate"
    ) %>%
    mutate(
      method = factor(method, levels = c("SAND", "TLGLM")),
      source = factor(source, levels = proxy_names[seq_len(K)])
    )
}

plot_selection_counts <- function(sel_df, outfile) {
  long_df <- selection_count_long(sel_df)

  p <- ggplot(long_df, aes(method, n_selected)) +
    geom_boxplot(width = 0.55, outlier.alpha = 0.35) +
    geom_jitter(
      width = 0.12,
      height = 0,
      alpha = 0.25,
      size = 0.8,
      color = "red"
    ) +
    labs(
      x = NULL,
      y = "Number of selected sources",
      title = "Source-set size: SAND vs TLGLM"
    ) +
    theme_classic(base_size = 11)

  ggsave(outfile, p, width = 4.2, height = 3.6)
  invisible(p)
}

plot_selection_scatter <- function(sel_df, outfile) {
  p <- ggplot(sel_df, aes(tlglm_n, sand_n)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(alpha = 0.45, size = 1.6, color = "red") +
    labs(
      x = "TLGLM selected sources",
      y = "SAND selected sources",
      title = "Paired source counts per CV split"
    ) +
    theme_classic(base_size = 11)

  ggsave(outfile, p, width = 4.2, height = 3.8)
  invisible(p)
}

plot_selection_overlap <- function(sel_df, outfile) {
  summary_df <- sel_df %>%
    summarise(
      jaccard_mean = mean(jaccard, na.rm = TRUE),
      jaccard_median = median(jaccard, na.rm = TRUE),
      exact_match_rate = mean(exact_match, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )

  p <- ggplot(sel_df, aes(jaccard)) +
    geom_histogram(binwidth = 0.05, fill = "grey75", color = "white") +
    geom_vline(
      xintercept = summary_df$jaccard_mean,
      linetype = 2,
      linewidth = 0.6
    ) +
    annotate(
      "text",
      x = summary_df$jaccard_mean,
      y = Inf,
      vjust = 1.4,
      hjust = -0.05,
      size = 3.2,
      label = sprintf(
        "mean = %.2f\nexact match = %.1f%%",
        summary_df$jaccard_mean,
        100 * summary_df$exact_match_rate
      )
    ) +
    labs(
      x = "Jaccard index",
      y = "Count",
      title = "Overlap between SAND and TLGLM selections"
    ) +
    theme_classic(base_size = 11)

  ggsave(outfile, p, width = 4.2, height = 3.6)
  invisible(list(plot = p, summary = summary_df))
}

plot_selection_frequency <- function(freq_df, outfile) {
  p <- ggplot(freq_df, aes(source, selection_rate, fill = method)) +
    geom_col(position = "dodge", width = 0.72) +
    scale_fill_manual(values = c(SAND = "grey35", TLGLM = "grey75")) +
    labs(
      x = NULL,
      y = "Selection frequency across CV splits",
      fill = NULL,
      title = "Per-source selection rate"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 8)
    )

  ggsave(outfile, p, width = 6.5, height = 3.8)
  invisible(p)
}

summarize_selection <- function(sel_df) {
  sel_df %>%
    summarise(
      n_splits = dplyr::n(),
      sand_n_mean = mean(sand_n, na.rm = TRUE),
      sand_n_sd = sd(sand_n, na.rm = TRUE),
      tlglm_n_mean = mean(tlglm_n, na.rm = TRUE),
      tlglm_n_sd = sd(tlglm_n, na.rm = TRUE),
      jaccard_mean = mean(jaccard, na.rm = TRUE),
      jaccard_sd = sd(jaccard, na.rm = TRUE),
      exact_match_rate = mean(exact_match, na.rm = TRUE),
      .groups = "drop"
    )
}

# --- run --------------------------------------------------------------------

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("cache", showWarnings = FALSE, recursive = TRUE)

wide_df <- load_pilot_results()
proxy_names <- infer_proxy_names()
K <- if (!is.null(proxy_names)) length(proxy_names) else NULL

sel_df <- enrich_selection(wide_df, K = K)
sel_summary <- summarize_selection(sel_df)
freq_df <- selection_frequency(sel_df, proxy_names = proxy_names)

saveRDS(wide_df, "cache/resdf_all.rds")
write.csv(wide_df, "cache/resdf_all.csv", row.names = FALSE)
write.csv(sel_summary, "cache/selection_summary.csv", row.names = FALSE)
write.csv(freq_df, "cache/selection_frequency.csv", row.names = FALSE)

plot_selection_counts(sel_df, "figures/source_selection_counts.pdf")
plot_selection_scatter(sel_df, "figures/source_selection_scatter.pdf")
overlap_out <- plot_selection_overlap(sel_df, "figures/source_selection_overlap.pdf")
plot_selection_frequency(freq_df, "figures/source_selection_frequency.pdf")

cat(
  "\nAggregated", nrow(wide_df), "CV rows from",
  length(unique(wide_df$source_file)), "cache file(s).\n\n"
)
cat("Selection summary:\n")
print(as.data.frame(sel_summary), row.names = FALSE)

cat("\nSaved:\n")
cat("  cache/resdf_all.rds\n")
cat("  cache/resdf_all.csv\n")
cat("  cache/selection_summary.csv\n")
cat("  cache/selection_frequency.csv\n")
cat("  figures/source_selection_counts.pdf\n")
cat("  figures/source_selection_scatter.pdf\n")
cat("  figures/source_selection_overlap.pdf\n")
cat("  figures/source_selection_frequency.pdf\n")
