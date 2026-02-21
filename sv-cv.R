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

sieve_solver <- function (model, Y, l1 = TRUE, family = "gaussian", lambda = NULL){
  if (l1 == FALSE) {
    Phi <- model$Phi
    beta_hat <- least_square_C(Phi, Y)
    model$beta_hat <- beta_hat
  } else {
    glmnet.fit <- glmnet(model$Phi[, -1], Y, family = family, alpha = 1, 
                         lambda = lambda, intercept = TRUE, standardize = FALSE)
    model$lambda <- glmnet.fit$lambda
    model$beta_hat <- coef(glmnet.fit)
    model$family <- family
  }
  return(model)
}

cross_validated_sieve <- function(X, Y, basis_numbers = NULL, n_folds = 5, type = "cosine") {
  p <- ncol(X)
  n_sample <- nrow(X)
  if (is.null(basis_numbers)) {
    basis_numbers = ceiling(c(
      p * c(5, n_sample^(1 / 5), n_sample^(1 / 3))
      , p^2 * c(5, n_sample^(1 / 5), n_sample^(1 / 3))
    ))
  }

  validation_split_idx = index_splitter(seq_len(n_sample), n_folds = n_folds)
  total_iterations = length(basis_numbers) * n_folds
  tuning_rows_list = vector("list", total_iterations)
  list_idx <- 1
  for (basis_n in basis_numbers) {
    global_model = sieve_preprocess(X, basis_n, type = type)
    global_fit = sieve_solver(global_model, Y, lambda = NULL)
    global_lambda = global_fit$lambda
    for (split_idx in seq_len(n_folds)) {
      train_X = X[-validation_split_idx[[split_idx]], , drop = FALSE]
      train_Y = Y[-validation_split_idx[[split_idx]]]
      validation_X = X[validation_split_idx[[split_idx]], , drop = FALSE]
      validation_Y = Y[validation_split_idx[[split_idx]]]
      sieve_model = sieve_preprocess(train_X, basisN = basis_n, type = type)
      sieve_fit = sieve_solver(sieve_model, train_Y, lambda = global_lambda)
      sieve_validation = sieve_predict(
        model = sieve_fit,
        testX = validation_X,
        testY = validation_Y
      )
      tuning_rows_list[[list_idx]] = data.table::data.table(
        l1_penalty = sieve_fit$lambda,
        l1_penalty_index = seq_along(sieve_fit$lambda),
        basis_n = rep(basis_n, length(sieve_fit$lambda)),
        validation_mse = sieve_validation$MSE,
        split_idx = rep(split_idx, length(sieve_fit$lambda))
      )
      list_idx <- list_idx + 1
    }
  }
  
  tuning_dt = data.table::rbindlist(tuning_rows_list)
  average_dt = tuning_dt[, .(average_mse = mean(validation_mse),
                             l1_penalty = mean(l1_penalty)), 
                         by = .(l1_penalty_index, basis_n)]

  best_idx = which.min(average_dt$average_mse)
  list(best_lambda = average_dt$l1_penalty[best_idx],
       best_basis_number = average_dt$basis_n[best_idx])
}

fit_sieve_with_cv <- function(X, Y, type = "cosine", n_folds = 5) {
  p = ncol(X)
  # basis_numbers <- p * c(2, 4, 8) # temporarily fixed for pilot study
  basis_numbers <- NULL
  best_hyperparameter = cross_validated_sieve(X, Y, basis_numbers, n_folds, type)
  sieve_model = sieve_preprocess(X, 
                                 basisN = best_hyperparameter$best_basis_number,
                                 type = type)
  sieve_fit = sieve_solver(sieve_model, Y, lambda = best_hyperparameter$best_lambda)
  list(hyperparameter = best_hyperparameter, sieve_fit = sieve_fit)
}

fmap_predict <- function(fmap, X, type = 'sieve') {
  n <- nrow(X)
  p2 <- length(fmap)
  hat_Z <- matrix(0, n, p2)
  if (type == 'sieve') {
    for (p2id in 1:p2) {
      hat_Z[, p2id] = sieve_predict(fmap[[p2id]]$sieve_fit, X)$predictY
    }
  } else if (type == 'linear') {
    for (p2id in 1:p2) {
      hat_Z[, p2id] = predict(fmap[[p2id]], X, s = 'lambda.min')
    }
  }
  colnames(hat_Z) <- names(fmap) 
  return(hat_Z)
}

cv_fitH_from_prxy <- function(x_common, x_prxy_only, type = "linear", n_folds = 5, n_cores = 1) {
  p2 <- ncol(x_prxy_only)
  fit_one_prxy_col <- function(p2id, type) {
    y <- x_prxy_only[, p2id]
    if (type == 'sieve') fit_sieve_with_cv(x_common, y, 'cosine', n_folds)
    else if (type == 'linear')
      glmnet(x_common, y, alpha = 0, lambda = 0)
      # cv.glmnet(x_common, y, alpha = 0, nfolds = n_folds)
  }
  if (.Platform$OS.type == "unix" && n_cores > 1) {
    if (requireNamespace("pbmcapply", quietly = TRUE)) {
      sieve_fits = pbmcapply::pbmclapply(
        seq_len(p2), fit_one_prxy_col, mc.cores = n_cores, type = type
      )
    } else {
      message("Install 'pbmcapply' to show a parallel progress bar: install.packages('pbmcapply')")
      sieve_fits = parallel::mclapply(seq_len(p2), fit_one_prxy_col, 
                                      type = type, mc.cores = n_cores)
    }
  } else sieve_fits = lapply(seq_len(p2), fit_one_prxy_col, type = type)
  max_len = max(lengths(sieve_fits))
  names(sieve_fits) <- colnames(x_prxy_only)
  sieve_fits
}