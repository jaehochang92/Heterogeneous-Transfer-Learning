## Modified from:
## https://github.com/terrytianyuzhang/Sieve/blob/master/R/cv_sieve_need_to_include.R
## Zhang, T., & Simon, N. (2023). Regression in tensor product spaces by the method of sieves. Electronic journal of statistics, 17(2), 3660.


library(data.table)
library(Sieve)

index_splitter <- function(index_array, n_folds = 5) {
  part_length = length(index_array) %/% n_folds
  fold_parts = vector("list", n_folds)
  shuffled_index = sample(index_array)
  for (fold_idx in seq_len(n_folds)) {
    start_idx = (fold_idx - 1) * part_length + 1
    end_idx = if (fold_idx < n_folds) fold_idx * part_length else length(index_array)
    fold_parts[[fold_idx]] = shuffled_index[start_idx:end_idx]
  }
  fold_parts
}

cross_validated_sieve <- function(all_data, basis_numbers = NULL, n_folds = 5, type = "cosine") {
  p <- ncol(all_data) - 1
  if (is.null(basis_numbers)) {
    basis_numbers = ceiling(c(
      p * c(5, nrow(all_data)^(1 / 5), nrow(all_data)^(1 / 3))
      # , p^2 * c(5, nrow(all_data)^(1 / 5), nrow(all_data)^(1 / 3))
    ))
  }

  validation_split_idx = index_splitter(seq_len(nrow(all_data)), n_folds = n_folds)
  tuning_rows = data.frame()

  for (basis_n in basis_numbers) {
    sieve_fitting_lambda = NULL
    for (split_idx in seq_len(n_folds)) {
      train_data = all_data[-validation_split_idx[[split_idx]], ]
      validation_data = all_data[validation_split_idx[[split_idx]], ]

      sieve_model = sieve_preprocess(X = train_data[, 2:(p + 1)], basisN = basis_n, type = type)
      sieve_fit = sieve_solver(model = sieve_model, Y = train_data$Y, lambda = sieve_fitting_lambda)
      sieve_fitting_lambda = sieve_fit$lambda

      sieve_validation = sieve_predict(
        model = sieve_fit,
        testX = validation_data[, 2:(p + 1)],
        testY = validation_data$Y
      )

      tuning_rows = rbind(
        tuning_rows,
        data.frame(
          l1_penalty = sieve_fit$lambda,
          l1_penalty_index = seq_along(sieve_fit$lambda),
          basis_n = rep(basis_n, length(sieve_fit$lambda)),
          validation_mse = sieve_validation$MSE,
          split_idx = rep(split_idx, length(sieve_fit$lambda))
        )
      )
    }
  }

  tuning_dt = data.table::data.table(tuning_rows)
  average_dt = tuning_dt[, .(
    average_mse = mean(validation_mse),
    l1_penalty = mean(l1_penalty)
  ), by = .(l1_penalty_index, basis_n)]

  best_idx = which.min(average_dt$average_mse)
  list(
    best_lambda = average_dt$l1_penalty[best_idx],
    best_basis_number = average_dt$basis_n[best_idx]
  )
}

fit_sieve_with_cv <- function(all_data, type = "cosine", n_folds = 5) {
  p = ncol(all_data) - 1 # first column is y
  basis_numbers <- NULL
  best_hyperparameter = cross_validated_sieve(
    all_data = all_data,
    basis_numbers = basis_numbers,
    n_folds = n_folds,
    type = type
  )
  sieve_model = sieve_preprocess(
    X = all_data[, 2:(p + 1)],
    basisN = best_hyperparameter$best_basis_number,
    type = type
  )
  sieve_fit = sieve_solver(
    model = sieve_model,
    Y = all_data$Y,
    lambda = best_hyperparameter$best_lambda
  )
  list(hyperparameter = best_hyperparameter, sieve_fit = sieve_fit)
}

# x_proxy: n x p matrix. First p1 features are common in proxy and target domains; remaining p2 are proxy-only.
cv_fitH_from_prxy <- function(x_proxy, p1, type = "linear", n_folds = 5, n_cores = 1) {
  p2 <- ncol(x_proxy) - p1
  if (p2 <= 0) stop("x_proxy must contain proxy-only columns beyond p1.")
  x_common = x_proxy[, 1:p1, drop = FALSE]
  x_proxy_only = x_proxy[, (p1 + 1):(p1 + p2), drop = FALSE]
  coef_list = vector("list", p2)
  if (type == 'sieve') {
    fit_one_proxy_col <- function(p2id) {
      df = data.frame(Y = x_proxy_only[, p2id], x_common, check.names = FALSE)
      sieve_cv_fit = fit_sieve_with_cv(df, type = 'cosine', n_folds = n_folds)
      as.numeric(sieve_cv_fit$sieve_fit$beta_hat)
    }
    if (.Platform$OS.type == "unix" && n_cores > 1) {
      if (requireNamespace("pbmcapply", quietly = TRUE)) {
        coef_list = pbmcapply::pbmclapply(
          seq_len(p2),
          fit_one_proxy_col,
          mc.cores = n_cores
        )
      } else {
        message("Install 'pbmcapply' to show a parallel progress bar: install.packages('pbmcapply')")
        coef_list = parallel::mclapply(seq_len(p2), fit_one_proxy_col, mc.cores = n_cores)
      }
    } else {
      coef_list = lapply(seq_len(p2), fit_one_proxy_col)
    }
    max_len = max(lengths(coef_list))
    emap = matrix(0, nrow = max_len, ncol = p2)
    for (p2id in seq_along(coef_list)) {
      coef_vec = coef_list[[p2id]]
      emap[seq_along(coef_vec), p2id] = coef_vec
    }
    colnames(emap) = colnames(x_proxy_only)
  } else if (type == 'linear') {
    lm_fit <- lm(x_proxy_only ~ x_common - 1)
    emap <- coef(lm_fit)
  }
  emap
}