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

setwd('~/git/Heterogeneous-Transfer-Learning-code/brewing../updated')
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
  make_option(c("--np_values"), type = "character", default = "200,400,800,1600"),
  make_option(c("--nt_values"), type = "character", default = "30,60,90"),
  make_option(c("--p1_values"), type = "character", default = "10, 20, 30"),
  make_option(c("--p2_values"), type = "character", default = "10, 20, 30")
)
opt <- parse_args(OptionParser(option_list = option_list))

slurm_array_task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
sweep_task_id <- as.integer(Sys.getenv("SIM_SWEEP_TASK_ID"))
if (is.na(sweep_task_id)) {
  sweep_task_id <- slurm_array_task_id
}
if (is.na(sweep_task_id)) {
  sweep_task_id <- 1L
}

global_task_id <- as.integer(Sys.getenv("SIM_GLOBAL_TASK_ID"))
if (is.na(global_task_id)) {
  global_task_id <- slurm_array_task_id
}
if (is.na(global_task_id)) {
  global_task_id <- sweep_task_id
}

log_dir <- "logs"
if (!dir.exists(log_dir))
  dir.create(log_dir)
log_path <- sprintf("%s/progress_task_%04d.log", log_dir, global_task_id)
log_message <- function(msg) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", stamp, msg),
      file = log_path,
      append = TRUE)
}

sweep_plan <- build_sweep_plan(opt)
total_tasks <- nrow(sweep_plan)
if (sweep_task_id > total_tasks) {
  log_message(sprintf(
    "Task id %d is out of range (max %d). Exiting cleanly.",
    sweep_task_id,
    total_tasks
  ))
  quit(save = "no", status = 0)
}

task_cfg <- sweep_plan[sweep_task_id]
sp_cases <- c("sp", "nsp")
fmap_cases <- c("linear", "nlinear")
total_steps <- opt$rep * length(sp_cases) * length(fmap_cases)
log_message(
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

result_list <- vector("list", total_steps)
row_idx <- 1L
for (i in seq_len(opt$rep)) {
  current_seed <- opt$seedno + i - 1L
  for (sp_case in sp_cases) {
    for (fmap in fmap_cases) {
      t0 <- proc.time()[3]
      log_message(
        sprintf(
          "Progress %d/%d: iter=%d sp_case=%s dgp_fmap=%s seed=%d",
          row_idx, total_steps, i, sp_case, fmap, current_seed
        )
      )
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
      if (row_idx %% 2L == 0L || row_idx == total_steps) {
        partial_dt <- rbindlist(result_list[seq_len(row_idx)], fill = TRUE)
        partial_path <- sprintf("results/sim_res_task_%04d.partial.rds", global_task_id)
        saveRDS(partial_dt, partial_path)
      }
      row_idx <- row_idx + 1L
    }
  }
}

final_dt <- rbindlist(result_list, fill = TRUE, ignore.attr = TRUE)

out_dir <- "results"
if (!dir.exists(out_dir))
  dir.create(out_dir)

save_path <- sprintf("%s/sim_res_task_%04d.rds", out_dir, global_task_id)
saveRDS(final_dt, save_path)
log_message(sprintf(
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