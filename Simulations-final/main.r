# main.R
library(data.table)
library(Matrix)
library(optparse)
library(MASS)
library(glmnet)
library(stargazer)
library(pbapply)
library(dplyr)
library(truncnorm)
library(Sieve)

# setwd('~/git/Heterogeneous-Transfer-Learning-code/Simulations-final')
source("src/methods.r")
source("src/dgp.r")
source("src/mapping.R")
source("src/htl.R")

# 1. simulation configurations (admits SLURM arguments)
option_list <- list(
  make_option(c("--seedno"), type = "integer", default = 1992),
  make_option(c("--K"), type = "integer", default = 1),
  make_option(c("--np"), type = "integer", default = 2e2),
  make_option(c("--nt"), type = "integer", default = 30),
  make_option(c("--p1"), type = "integer", default = 20),
  make_option(c("--p2"), type = "integer", default = 20),
  make_option(c("--rep"), type = "integer", default = 10),
  make_option(c("--K_values"), type = "character", default = "2,4,8,16"),
  make_option(c("--np_values"), type = "character", default = "200,400,800"),
  make_option(c("--nt_values"), type = "character", default = "30,60,90"),
  make_option(c("--p1_values"), type = "character", default = "10, 20, 30"),
  make_option(c("--p2_values"), type = "character", default = "10, 20, 30"),
  make_option(c("--sweep_params"), type = "character", default = "K,np,nt,p1,p2"),
  make_option(c("--out_dir"), type = "character", default = "results"),
  make_option(c("--log_dir"), type = "character", default = "rlogs"),
  make_option(c("--write_progress_log"), action = "store_true", default = FALSE),
  make_option(c("--verbose"), action = "store_true", default = FALSE),
  make_option(c("--save_partial"), action = "store_true", default = FALSE),
  make_option(c("--partial_every"), type = "integer", default = 0),
  make_option(c("--keep_partial"), action = "store_true", default = FALSE),
  make_option(c("--bundle_id"), type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

global_task_id <- as.integer(Sys.getenv("SIM_GLOBAL_TASK_ID"))
slurm_array_task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(global_task_id)) {
  global_task_id <- slurm_array_task_id
}
if (is.na(global_task_id)) {
  global_task_id <- 1L
}

bundle_id <- trimws(opt$bundle_id)
if (!nzchar(bundle_id)) {
  bundle_id <- trimws(Sys.getenv("SIM_BUNDLE_ID", unset = ""))
}
if (nzchar(bundle_id)) {
  bundle_id <- gsub("[^A-Za-z0-9._-]+", "-", bundle_id)
}

verbose <- isTRUE(opt$verbose)
write_progress_log <- isTRUE(opt$write_progress_log)

out_dir_base <- opt$out_dir
log_dir_base <- opt$log_dir
log_path <- NULL

log_info <- function(msg) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s\n", stamp, msg)
  cat(line)
  if (!is.null(log_path)) {
    cat(line, file = log_path, append = TRUE)
  }
}

log_progress <- function(msg) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s\n", stamp, msg)
  if (verbose) {
    cat(line)
  }
  if (!is.null(log_path)) {
    cat(line, file = log_path, append = TRUE)
  }
}

sweep_plan <- build_sweep_plan(opt)
sweep_params <- trimws(unlist(strsplit(opt$sweep_params, ",", fixed = TRUE)))
sweep_params <- sweep_params[nzchar(sweep_params)]
if (length(sweep_params) == 0L) {
  stop("No sweep parameters selected. Provide --sweep_params with at least one of K,np,nt,p1,p2")
}
sweep_plan <- sweep_plan[sweep_param %in% sweep_params]
sweep_plan <- unique(sweep_plan)
total_tasks <- nrow(sweep_plan)
sweep_task_id <- ((global_task_id - 1L) %% total_tasks) + 1L
if (sweep_task_id > total_tasks) {
  log_info(sprintf(
    "Task id %d is out of range (max %d). Exiting cleanly.",
    sweep_task_id,
    total_tasks
  ))
  quit(save = "no", status = 0)
}

task_cfg <- sweep_plan[sweep_task_id]
sweep_param_dir <- gsub("[^A-Za-z0-9._-]+", "-", as.character(task_cfg$sweep_param))
if (!nzchar(sweep_param_dir)) {
  sweep_param_dir <- "unknown"
}

out_dir <- file.path(out_dir_base, sweep_param_dir)
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

if (write_progress_log) {
  log_dir <- file.path(log_dir_base, sweep_param_dir)
  if (!dir.exists(log_dir))
    dir.create(log_dir, recursive = TRUE)
  log_path <- if (nzchar(bundle_id)) {
    sprintf("%s/progress_bundle_%s.log", log_dir, bundle_id)
  } else {
    sprintf("%s/progress_task_%04d.log", log_dir, global_task_id)
  }
}

sp_cases <- c("sp", "nsp")
fmap_cases <- c("linear", "nlinear")
total_steps <- opt$rep * length(sp_cases) * length(fmap_cases)
log_info(
  sprintf(
    "Start global_task=%d sweep_task=%d/%d sweep_param=%s K=%d np=%d nt=%d p1=%d p2=%d rep=%d",
    global_task_id,
    sweep_task_id,
    total_tasks,
    task_cfg$sweep_param,
    task_cfg$K,
    task_cfg$np,
    task_cfg$nt,
    task_cfg$p1,
    task_cfg$p2,
    opt$rep
  )
)

save_partial <- isTRUE(opt$save_partial)
keep_partial <- isTRUE(opt$keep_partial)
partial_every <- as.integer(opt$partial_every)
if (save_partial && (is.na(partial_every) || partial_every <= 0L)) {
  partial_every <- 2L
}

save_path <- if (nzchar(bundle_id)) {
  sprintf("%s/sim_res_bundle_%s.rds", out_dir, bundle_id)
} else {
  sprintf("%s/sim_res_task_%04d.rds", out_dir, global_task_id)
}
partial_path <- if (save_partial) {
  if (nzchar(bundle_id)) {
    sprintf("%s/sim_res_bundle_%s.partial.rds", out_dir, bundle_id)
  } else {
    sprintf("%s/sim_res_task_%04d.partial.rds", out_dir, global_task_id)
  }
} else {
  NULL
}

result_list <- vector("list", total_steps)
row_idx <- 1L
for (i in seq_len(opt$rep)) {
  current_seed <- opt$seedno + i - 1L
  for (sp_case in sp_cases) {
    for (fmap in fmap_cases) {
      t0 <- proc.time()[3]
      ran_once <- sim.runner(
        seedno = current_seed,
        K = task_cfg$K,
        np = task_cfg$np,
        nt = task_cfg$nt,
        ntest = 100,
        p1 = task_cfg$p1,
        p2 = task_cfg$p2,
        dgp_fmap = fmap,
        sp = sp_case
      )
      elapsed_sec <- proc.time()[3] - t0
      est_err_named <- setNames(ran_once$Est_Error, paste0("Est_", names(ran_once$Est_Error)))
      pe_named <- setNames(ran_once$Prediction_Error, paste0("PE_", names(ran_once$Prediction_Error)))
      result_list[[row_idx]] <- c(
        list(
          Task = global_task_id,
          SweepTask = sweep_task_id,
          TotalTasks = total_tasks,
          SweepParam = task_cfg$sweep_param,
          K = task_cfg$K,
          np = task_cfg$np,
          nt = task_cfg$nt,
          p1 = task_cfg$p1,
          p2 = task_cfg$p2,
          dgp_fmap = fmap,
          Sp = sp_case,
          Iter = i,
          Seed = current_seed,
          RuntimeSec = elapsed_sec
        ),
        as.list(est_err_named),
        as.list(pe_named)
      )
      log_progress(
        sprintf(
          "Progress %d/%d: iter=%d sp_case=%s dgp_fmap=%s seed=%d",
          row_idx, total_steps, i, sp_case, fmap, current_seed
        )
      )
      if (save_partial && (row_idx %% partial_every == 0L || row_idx == total_steps)) {
        partial_dt <- rbindlist(result_list[seq_len(row_idx)], fill = TRUE, ignore.attr = TRUE)
        if (nzchar(bundle_id) && file.exists(save_path)) {
          existing_dt <- readRDS(save_path)
          partial_dt <- rbindlist(list(existing_dt, partial_dt), fill = TRUE, use.names = TRUE)
        }
        saveRDS(partial_dt, partial_path)
      }
      row_idx <- row_idx + 1L
    }
  }
}

final_dt <- rbindlist(result_list, fill = TRUE, ignore.attr = TRUE)

if (nzchar(bundle_id) && file.exists(save_path)) {
  existing_dt <- readRDS(save_path)
  final_dt <- rbindlist(list(existing_dt, final_dt), fill = TRUE, use.names = TRUE)
}

saveRDS(final_dt, save_path)
if (!is.null(partial_path) && file.exists(partial_path) && !keep_partial) {
  file.remove(partial_path)
}

log_info(sprintf(
  "Completed global_task=%d sweep_task=%d rows=%d saved=%s",
  global_task_id,
  sweep_task_id,
  nrow(final_dt),
  save_path
))
cat(sprintf(
  "Global task %d (sweep task %d) completed. Saved to %s\n",
  global_task_id,
  sweep_task_id,
  save_path
))