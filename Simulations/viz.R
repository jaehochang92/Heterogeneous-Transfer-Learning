options(expressions = 5e5)
suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(purrr)
  library(rlang)
  library(gridExtra)
})

source('codes/simulation.R')
# dev.off()

# ---------- helpers ----------------------------------------------------
safe_read_rds <- function(path) {
  tryCatch(readRDS(path), error = function(e) {
    warning(glue("Missing or unreadable: {path}"))
    NULL
  })
}

load_param <- function(param, values, dir = "rds") {
  p <- sym(param)
  map_dfr(values, function(v) {
    f <- glue("{dir}/{param}={v}.rds")
    obj <- safe_read_rds(f)
    if (is.null(obj)) return(NULL)
    bind_cols(obj, !!p := v)
  })
}

maybe_drop_lasso <- function(df, draw_lasso) {
  if (isTRUE(draw_lasso)) df else 
    dplyr::filter(df, .data$Methodology != "lasso_trgt")
}

make_two_plots <- function(df, param, levels_vec, title_str, out_prefix,
                           width, height, eng, quiet = FALSE) {
  # Empty or NULL df? — bail out quietly
  if (is.null(df) || nrow(df) == 0) {
    if (!quiet) message(glue("[{out_prefix}] skipped: empty df"))
    return(invisible(list(EstErr = NULL, PrdErr = NULL)))
  }
  
  # Param column must exist
  if (!param %in% names(df)) {
    if (!quiet) message(glue("[{out_prefix}] skipped: missing column `{param}`"))
    return(invisible(list(EstErr = NULL, PrdErr = NULL)))
  }
  
  p <- sym(param)
  
  base <- df %>%
    mutate(!!p := factor(.data[[param]], levels = levels_vec))
  
  build_plot <- function(dat, loss_name, ycol, ycalc) {
    bxplt <- geom_boxplot(aes(fill = Methodology), outliers = F, 
                          size = .3, outlier.size = 0.2)
    fctgrd <- facet_wrap(. ~ dlt_sprs, labeller = label_both, 
                         scales = 'free_x')
    scl_thm <- scale_fill_manual(values = c(
      "TL_impute" = "#FFD700",
      "TL_matched" = "#9ACD32",
      "lasso_trgt" = "#d1fEF3",
      "oracle" = "#FFB6C1"
    ))
    g_theme <- theme_bw() + theme(
      legend.position = "top",
      legend.justification = "right",
      legend.margin = margin(1, 1, 1, 1),
      legend.box.margin = margin(0, 0, -5, 0), 
      plot.margin = unit(rep(.3, 4), "cm"),
      strip.background = element_blank(),
      panel.spacing = unit(.5, "lines")
      # ,axis.title.y = element_blank()
    )
    
    dsub <- dat %>% filter(.data$loss == loss_name)
    if (nrow(dsub) == 0) {
      if (!quiet) message(glue("[{out_prefix}] no rows for loss = '{loss_name}'"))
      return(NULL)
    }
    gg <- ggplot(dsub, aes(x = !!p, y = value)) +
      bxplt + fctgrd + scl_thm + g_theme +
      ggtitle(title_str) + ylab(dsub$loss[1])
    gg
  }
  
  ggE <- build_plot(base, "EstErr", "EstErr", value)
  if (!is.null(ggE)) {
    ggsave(glue("{eng}/EstErr-{out_prefix}.{eng}"), ggE,
           width = width, height = height, units = "in")
  }
  
  ggP <- build_plot(base, "PrdErr", "PrdErr", value)
  if (!is.null(ggP)) {
    ggsave(glue("{eng}/PrdErr-{out_prefix}.{eng}"), ggP,
           width = width, height = height, units = "in")
  }
  gga <- grid.arrange(ggP, ggE, nrow = 1)
  print(gga)
  invisible(list(EstErr = ggE, PrdErr = ggP))
}

# Config ------------------------------------------------------------------

eng <- 'eps'
eng <- 'png'

# ---------- global output dir -----------------------------------------
dir.create(eng, showWarnings = FALSE, recursive = TRUE)

fcol <- 2
frow <- 1
wd <- ht <- 2.5 # plot width and height
bed <- .75 # bedding for height

.plot_w <- wd * fcol
.plot_h <- bed + ht * frow

draw_lasso <- T


dfK <- load_param("K", Kv) %>% maybe_drop_lasso(draw_lasso)
title_K <- glue("(np, nt, ntest, p1, p2) = ({np0}, {nt0}, {ntest}, {p10}, {p20})")
make_two_plots(dfK, "K", Kv, title_K, "K", .plot_w, .plot_h, eng)

dfnt <- load_param("nt", ntv) %>% maybe_drop_lasso(draw_lasso)
title_nt <- glue("(K, np, ntest, p1, p2) = ({1}, {np0}, {ntest}, {p10}, {p20})")
make_two_plots(dfnt, "nt", ntv, title_nt, "nt", .plot_w, .plot_h, eng)

dfnp <- load_param("np", npv) %>% maybe_drop_lasso(draw_lasso)
title_np <- glue("(K, nt, ntest, p1, p2) = ({1}, {nt0}, {ntest}, {p10}, {p20})")
make_two_plots(dfnp, "np", npv, title_np, "np", .plot_w, .plot_h, eng)

dfp1 <- load_param("p1", p1v) %>% maybe_drop_lasso(draw_lasso)
title_p1 <- glue("(K, np, nt, ntest, p2) = ({1}, {np0}, {nt0}, {ntest}, {p20})")
make_two_plots(dfp1, "p1", p1v, title_p1, "p1", .plot_w, .plot_h, eng)

dfp2 <- load_param("p2", p2v) %>% maybe_drop_lasso(draw_lasso)
title_p2 <- glue("(K, np, nt, ntest, p1) = ({1}, {np0}, {nt0}, {ntest}, {p10})")
make_two_plots(dfp2, "p2", p2v, title_p2, "p2", .plot_w, .plot_h, eng)