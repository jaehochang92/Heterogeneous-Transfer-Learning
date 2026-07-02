library(data.table)
library(Matrix)
library(MASS)
library(glmnet)
library(pbapply)
library(dplyr)
library(glmtrans)
library(scalreg)

pilot_dir <- Sys.getenv(
  "SLURM_SUBMIT_DIR",
  unset = normalizePath(".", winslash = "/")
)
setwd(pilot_dir)

source("methods/methods.r")
source("methods/mapping.R")
source("methods/htl.R")
source("methods/sand.R")

assign_replicate_range <- function(array_task_id, submitted_tasks, n_replicates) {
  base_load <- n_replicates %/% submitted_tasks
  extra <- n_replicates %% submitted_tasks

  if (array_task_id <= extra) {
    my_load <- base_load + 1L
    start_idx <- (array_task_id - 1L) * (base_load + 1L) + 1L
  } else {
    my_load <- base_load
    start_idx <- extra * (base_load + 1L) + (array_task_id - extra - 1L) * base_load + 1L
  }

  if (my_load <= 0L) {
    return(integer())
  }
  seq.int(start_idx, start_idx + my_load - 1L)
}

array_task_id <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
submitted_tasks <- suppressWarnings(as.integer(Sys.getenv("SUBMITTED_TASKS", "1")))
n_replicates <- suppressWarnings(as.integer(Sys.getenv("N_REPLICATES", "1")))
seed_base <- suppressWarnings(as.integer(Sys.getenv("SEED_BASE", "2175")))
n_folds <- suppressWarnings(as.integer(Sys.getenv("N_FOLDS", "5")))
bundle_id <- Sys.getenv("SIM_BUNDLE_ID", paste0("local_", array_task_id))

if (any(is.na(c(array_task_id, submitted_tasks, n_replicates, seed_base, n_folds)))) {
  stop("Invalid SLURM / replicate environment variables.", call. = FALSE)
}

replicate_ids <- assign_replicate_range(array_task_id, submitted_tasks, n_replicates)
if (!length(replicate_ids)) {
  message("Task ", array_task_id, ": no assigned replicates; exiting.")
  quit(save = "no", status = 0)
}

target <- readRDS("target.rds")
proxy_list <- readRDS("proxy.rds")

slurm_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(slurm_cores) || slurm_cores < 1L) {
  slurm_cores <- 1L
}

res_df <- list()

for (rep_id in replicate_ids) {
  set.seed(seed_base + rep_id)
  fold_ids <- .index_spliter(seq_along(target$y), n_folds)

  for (fold_idx in seq_along(fold_ids)) {
    tst_id <- fold_ids[[fold_idx]]
    train <- test <- target
    train$X <- train$X[-tst_id, , drop = FALSE]
    train$y <- train$y[-tst_id]
    test$X <- test$X[tst_id, , drop = FALSE]
    test$y <- test$y[tst_id]

    train.tlglm <- list(x = as.matrix(train$X), y = as.numeric(train$y))
    proxy_list.tlglm <- lapply(proxy_list, function(p) {
      list(x = as.matrix(p$D$X), y = as.numeric(p$y))
    })

    tlglm_fit <- glmtrans(
      target = train.tlglm,
      source = proxy_list.tlglm,
      family = "gaussian",
      transfer.source.id = "auto",
      nfolds = 5,
      cores = slurm_cores
    )
    h_ytlglm_test <- predict(tlglm_fit, newx = test$X)
    pe_tlglm.org <- .calc_rmse(h_ytlglm_test, test$y)

    proxy_fit_all <- agg.proxy(proxy_list, "linear")
    res_all_htl <- htl(proxy_fit_all, train, test)

    proxy_preds <- make_proxy_preds_linear(train, proxy_list)
    hA <- run_sand(
      Xt = train$X,
      yt = train$y,
      proxy_preds_list = proxy_preds,
      gamma_scale = 1
    )

    if (length(hA) > 0L) {
      proxy_list_selected <- proxy_list[hA]
      proxy_fit_sand_htl <- agg.proxy(proxy_list_selected, fmap = "linear")
      res_sand_htl <- htl(proxy_fit_sand_htl, train, test)
    } else {
      base_fit <- glmtrans(
        target = train.tlglm,
        family = "gaussian",
        transfer.source.id = 0,
        nfolds = 5,
        cores = slurm_cores
      )
      h_y_base_test <- predict(base_fit, newx = test$X)
      pe_base <- .calc_rmse(h_y_base_test, test$y)
      res_sand_htl <- list(pe = pe_base)
    }

    lasso_fit <- cv.glmnet(train$X, train$y, alpha = 1)
    lasso_pred <- predict(lasso_fit, newx = test$X, s = "lambda.min")
    lasso_pe <- .calc_rmse(test$y, lasso_pred)

    sand_ids <- sort(as.integer(hA))
    tlglm_ids <- tlglm_fit$transfer.source.id
    if (identical(tlglm_ids, "all")) {
      tlglm_ids <- seq_along(proxy_list)
    }
    tlglm_ids <- sort(as.integer(tlglm_ids))

    res_df[[length(res_df) + 1L]] <- data.frame(
      replicate = rep_id,
      fold = fold_idx,
      ALL_HTL = res_all_htl$pe,
      TLGLM_HmTL = pe_tlglm.org,
      SAND_HTL = res_sand_htl$pe,
      Lasso = lasso_pe,
      sand_sel = I(list(sand_ids)),
      tlglm_sel = I(list(tlglm_ids)),
      sand_n = length(sand_ids),
      tlglm_n = length(tlglm_ids),
      stringsAsFactors = FALSE
    )
  }
}

res_df <- bind_rows(res_df)

dir.create("cache", showWarnings = FALSE, recursive = TRUE)
out_path <- file.path("cache", paste0("resdf_", bundle_id, ".rds"))
saveRDS(res_df, out_path)

cat(
  "Task", array_task_id, ":",
  length(replicate_ids), "replicate(s),",
  nrow(res_df), "CV rows saved to", out_path, "\n"
)
