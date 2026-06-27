# main_sand.R
# =============================================================================
# SAND (Signal-Based Adaptive Negative Transfer Defense) Simulation Runner
# Linear feature mapping case only
#
# Usage (interactive):
#   Rscript main_sand.R
#
# Usage (SLURM / CLI with arguments):
#   Rscript main_sand.R --K 4 --np 400 --nt 50 --p1 20 --p2 20 --rep 50
#
# This script evaluates SAND source selection accuracy alongside HTL
# estimation / prediction errors across varying configurations.
# =============================================================================

library(data.table)
library(Matrix)
library(optparse)
library(MASS)
library(glmnet)
library(pbapply)
library(dplyr)
library(glmtrans)
library(scalreg)
library(xgboost)

source("src/methods.r")
source("src/dgp.r")
source("src/mapping.R")
source("src/htl.R")
source("src/sand.R")

# ---------------------------------------------------------------------------
# 1. CLI / SLURM argument parsing
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--seedno"),      type = "integer",   default = 1992),
  make_option(c("--h"),           type = "double",   default = 0.3,
              help = "Delta_P [default %default]"),
  make_option(c("--K"),           type = "integer",   default = 4,
              help = "Number of proxy sources [default %default]"),
  make_option(c("--np"),          type = "integer",   default = 100),
  make_option(c("--nt"),          type = "integer",   default = 30),
  make_option(c("--p1"),          type = "integer",   default = 20),
  make_option(c("--p2"),          type = "integer",   default = 20),
  make_option(c("--rep"),         type = "integer",   default = 20,
              help = "Number of Monte Carlo replicates [default %default]"),
  make_option(c("--h_values"),    type = "character", default = "0.3,0.6,1.2,2.4"),
  make_option(c("--K_values"),    type = "character", default = "2,4,8"), 
  # K_values should be even
  make_option(c("--np_values"),   type = "character", default = "200,400,800"),
  make_option(c("--nt_values"),   type = "character", default = "30,60,90"),
  make_option(c("--p1_values"),   type = "character", default = "10,20,30"),
  make_option(c("--p2_values"),   type = "character", default = "10,20,30"),
  make_option(c("--sweep_params"),type = "character", default = "h,K,np,nt,p1,p2",
              help = "Comma-separated sweep axes [default %default]"),
  make_option(c("--gamma_scale"), type = "double",    default = 1.0,
              help = "Multiplicative scale for gamma_n [default %default]"),
  make_option(c("--out_dir"),     type = "character", default = "results_sand"),
  make_option(c("--log_dir"),     type = "character", default = "rlogs_sand"),
  make_option(c("--verbose"),     action = "store_true", default = FALSE),
  make_option(c("--bundle_id"),   type = "character", default = ""),
  make_option(c("--task_start"),  type = "integer", default = 1L),
  make_option(c("--task_count"),  type = "integer", default = 1L)
)
opt <- parse_args(OptionParser(option_list = option_list))

array_task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(array_task_id)) array_task_id <- 1L

bundle_id <- trimws(opt$bundle_id)
if (!nzchar(bundle_id)) bundle_id <- trimws(Sys.getenv("SIM_BUNDLE_ID", unset = ""))
if (nzchar(bundle_id)) bundle_id <- gsub("[^A-Za-z0-9._-]+", "-", bundle_id)

verbose <- isTRUE(opt$verbose)

log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}
log_progress <- function(msg) {
  if (verbose) log_info(msg)
}

# ---------------------------------------------------------------------------
# 2. Build sweep plan (same helper as main.r)
# ---------------------------------------------------------------------------
sweep_plan   <- build_sweep_plan(opt)
sweep_params <- trimws(unlist(strsplit(opt$sweep_params, ",", fixed = TRUE)))
sweep_params <- sweep_params[nzchar(sweep_params)]
if (length(sweep_params) == 0L) stop("--sweep_params must contain at least one axis.")

sweep_plan   <- sweep_plan[sweep_param %in% sweep_params]
sweep_plan   <- unique(sweep_plan)
total_tasks  <- nrow(sweep_plan)

# ---------------------------------------------------------------------------
# 3. Single-replicate SAND runner
# ---------------------------------------------------------------------------
# sand_runner_once(1, 4, 2, 2e2, 30, 10, 10, 'linear')
sand_runner_once <- function(seedno, K, K_inf, np, nt, p1, p2, dgp_fmap,
                             sp = "nsp", gamma_scale = 1.0, h) {
  cnfg   <- dgp.confg(seedno, K, p1, p2, fmap = dgp_fmap, K_inf = K_inf, h = h)
  D_data <- dgp.gen.D(seedno, dgp_fmap, np, nt, ntest = 100, cnfg$Prxy, cnfg$Trgt)
  
  target <- D_data$targetdata
  test   <- D_data$testdata
  beta_t <- D_data$βt
  
  # --- Ground truth: inject adversarial sources and record true A --------
  true_A           <- cnfg$Prxy$true_A 
  proxy_list <- D_data$proxydata[[sp]]
  
  # --- SAND: source selection --------------------------------------------
  proxy_preds <- make_proxy_preds_linear(target, proxy_list)
  
  hA <- run_sand(
    X_t              = target$X,
    y_t              = target$Y,
    proxy_preds_list = proxy_preds,
    gamma_scale      = gamma_scale
  )
  
  # --- Selection accuracy metrics ----------------------------------------
  K_total       <- K
  n_selected    <- length(hA)
  true_pos      <- length(intersect(hA, true_A))
  false_pos     <- length(setdiff(hA, true_A))
  false_neg     <- length(setdiff(true_A, hA))
  precision     <- if (n_selected > 0L) true_pos / n_selected else 0.0
  recall        <- if (length(true_A) > 0L) true_pos / length(true_A) else 0.0
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0.0
  exact_match   <- as.integer(setequal(hA, true_A))
  
  # --- Baseline: glmtrans framework (Tian & Feng, 2023)
  target.tlglm <- list(x = as.matrix(target$X), y = as.numeric(target$Y))
  proxy_list.tlglm <- lapply(proxy_list, function(p) list(x = as.matrix(p$D$X), y = as.numeric(p$Y)))
  slurm_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  if (is.na(slurm_cores)) slurm_cores <- 1L
  tlglm_fit <- glmtrans(
    target             = target.tlglm,
    source             = proxy_list.tlglm,
    family             = "gaussian",
    transfer.source.id = "auto",
    nfolds             = 5,
    cores              = slurm_cores
  )
  h_y_tlglm_test    <- predict(tlglm_fit, newx = as.matrix(test$X))
  tlglm_coef        <- tlglm_fit$beta
  est_err_tlglm.org <- .calc_rmse(tlglm_coef[-1], beta_t[1:ncol(target$X)])
  pe_tlglm.org      <- .calc_rmse(h_y_tlglm_test, test$Y)
  
  # --- HTL with SAND-selected sources (if any selected) ------------------
  if (n_selected > 0L) {
    proxy_list_selected <- proxy_list[hA]
    proxy_fit_sand_htl <- agg.proxy(proxy_list_selected, fmap = "linear")
    res_sand_htl   <- htl.eval(proxy_fit_sand_htl, target, test, beta_t)
  } else {
    base_fit <- glmtrans(
    target             = target.tlglm,
    family             = "gaussian",
    transfer.source.id = 0,
    nfolds             = 5,
    cores              = slurm_cores
  )
  h_y_base_test    <- predict(base_fit, newx = as.matrix(test$X))
  base_coef        <- base_fit$beta
  est_err_base <- .calc_rmse(base_coef[-1], beta_t[1:ncol(target$X)])
  pe_base      <- .calc_rmse(h_y_base_test, test$Y)
    res_sand_htl <- list(est_err = est_err_base, pe = pe_base)
  }
  
  # --- Baseline: HTL using ALL sources -----------------------------------
  proxy_fit_all <- agg.proxy(proxy_list, "linear")
  res_all_htl   <- htl.eval(proxy_fit_all, target, test, beta_t)
  
  # --- Baseline: HmTL & HTL using only TRUE informative sources ----------
  proxy_fit_oracle <- agg.proxy(proxy_list[true_A], "none")
  res_oracle_hmtl  <- htl.eval(proxy_fit_oracle, target, test, beta_t)
  proxy_fit_oracle <- agg.proxy(proxy_list[true_A], "linear")
  res_oracle_htl   <- htl.eval(proxy_fit_oracle, target, test, beta_t)
  
  # --- Target-only baselines ---------------------------------------------
  res_lasso <- vanilla.glmnet(target, test, beta_t, 1)
  res_ridge <- vanilla.glmnet(target, test, beta_t, 0)
  
  result <- list(
    # Selection metrics
    Precision    = precision,
    Recall       = recall,
    F1           = f1,
    ExactMatch   = exact_match,
    #N_selected   = n_selected,
    #FalsePos     = false_pos,
    #FalseNeg     = false_neg,
    
    # Estimation errors
    # Est_Ridge      = res_ridge$est_err,
    # Est_Lasso      = res_lasso$est_err,
    # Est_AllHmTL   = res_all_hmtl$est_err,
    Est_TLGLM_org  = est_err_tlglm.org,
    Est_OracleHmTL = res_oracle_hmtl$est_err,
    Est_AllHTL    = res_all_htl$est_err,
    Est_SAND_HTL   = res_sand_htl$est_err,
    Est_OracleHTL  = res_oracle_htl$est_err,
    
    # Prediction errors
    # PE_Ridge       = res_ridge$pe,
    # PE_Lasso       = res_lasso$pe,
    # PE_AllHmTL    = res_all_hmtl$pe,
    PE_TLGLM_org   = pe_tlglm.org,
    PE_OracleHmTL  = res_oracle_hmtl$pe,
    PE_AllHTL     = res_all_htl$pe,
    PE_SAND_HTL    = res_sand_htl$pe,
    PE_OracleHTL   = res_oracle_htl$pe
  )
  result
}

# ---------------------------------------------------------------------------
# 4. Monte Carlo loop (Batch Processing)
# ---------------------------------------------------------------------------
fmap <- c("linear")
sp_cases <- c("sp", "nsp")

# Bash에서 전달받은 작업 범위 할당
task_start <- opt$task_start
task_count <- opt$task_count

all_results <- list()
row_idx <- 1L

for (offset in 0:(task_count - 1L)) {
  current_global_id <- task_start + offset
  sweep_task_id <- ((current_global_id - 1L) %% total_tasks) + 1L
  task_cfg      <- sweep_plan[sweep_task_id]
  repeat_idx    <- ((current_global_id - 1L) %/% total_tasks) + 1L
  current_seed  <- opt$seedno + repeat_idx - 1L
  log_progress(sprintf(
    ">>> Batch [%d/%d] | global=%d sweep=%d/%d param=%s h=%d K=%d np=%d nt=%d seed=%d",
    offset + 1L, task_count, current_global_id, sweep_task_id, total_tasks,
    task_cfg$sweep_param, task_cfg$h, task_cfg$K, task_cfg$np, task_cfg$nt, current_seed
  ))
  
  for (sp_case in sp_cases) {
    t0 <- proc.time()[3]
    res <- tryCatch(
      sand_runner_once(
        seedno      = current_seed,
        K           = task_cfg$K,
        K_inf       = as.integer((task_cfg$K)^(2/3)),
        np          = task_cfg$np,
        nt          = task_cfg$nt,
        p1          = task_cfg$p1,
        p2          = task_cfg$p2,
        dgp_fmap    = fmap,
        sp          = sp_case,
        gamma_scale = opt$gamma_scale,
        h           = task_cfg$h
      ),
      error = function(e) {
        log_info(sprintf("ERROR global=%d sp=%s seed=%d: %s", current_global_id, sp_case, current_seed, e$message))
        NULL
      }
    )
    elapsed <- proc.time()[3] - t0
    if (!is.null(res)) {
      all_results[[row_idx]] <- c(
        list(
          Task        = current_global_id,
          SweepTask   = sweep_task_id,
          TotalTasks  = total_tasks,
          SweepParam  = task_cfg$sweep_param,
          dgp_fmap    = fmap,
          h           = task_cfg$h,
          K           = task_cfg$K,
          K_inf       = as.integer((task_cfg$K)^(2/3)),
          np          = task_cfg$np,
          nt          = task_cfg$nt,
          p1          = task_cfg$p1,
          p2          = task_cfg$p2,
          Sp          = sp_case,
          Iter        = repeat_idx,
          Seed        = current_seed,
          RuntimeSec  = elapsed,
          gamma_scale = opt$gamma_scale
        ),
        res
      )
      row_idx <- row_idx + 1L
    }
  }
}

# ---------------------------------------------------------------------------
# 5. Save results
# ---------------------------------------------------------------------------
valid_rows <- Filter(Negate(is.null), all_results)
if (length(valid_rows) == 0L) {
  log_info("ERROR: No valid results were generated in this batch.")
  quit(save = "no", status = 1)
}

results_dt <- rbindlist(valid_rows, fill = TRUE, use.names = TRUE)

out_dir <- file.path(opt$out_dir, "batch_bundles")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save_path <- if (nzchar(bundle_id)) {
  file.path(out_dir, sprintf("sand_res_bundle_%s.rds", bundle_id))
} else {
  file.path(out_dir, sprintf("sand_res_batch_%04d_to_%04d.rds", task_start, task_start + task_count - 1L))
}

if (file.exists(save_path)) {
  existing_dt <- readRDS(save_path)
  results_dt <- rbindlist(list(existing_dt, results_dt), fill = TRUE, use.names = TRUE)
}

saveRDS(results_dt, save_path)
log_info(sprintf("Saved %d rows → %s", nrow(results_dt), save_path))

# ---------------------------------------------------------------------------
# 6. Print summary (interactive use)
# ---------------------------------------------------------------------------
if (verbose || array_task_id == 1L) {
  cols_sel  <- c("Precision", "Recall", "F1", "ExactMatch", "N_selected")
  cols_est  <- c()
  cols_pe   <- c("PE_SAND_HTL", 'PE_TLGLM_org')
  num_cols  <- intersect(c(cols_sel, cols_est, cols_pe), names(results_dt))
  
  summary_dt <- results_dt[,
                           lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 4)),
                           .SDcols = num_cols,
                           by = c("Sp", "SweepParam")
  ]
  cat("\n=== SAND Simulation Summary ===\n")
  print(summary_dt)
}