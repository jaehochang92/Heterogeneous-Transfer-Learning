library(data.table)
library(Sieve)
library(Matrix)  # Ensure the Matrix package is loaded for sparse matrices

# primitives --------------------------------------------------------------


calc_rmse = function(y_trth, y_prd) sqrt(mean((y_trth - y_prd)^2))

# Function to zero-pad and cbind two sparse matrices
padder0 <- function(...) {
  mats <- list(...)
  if (length(mats) == 0) return(NULL)
  mats <- lapply(mats, function(x) {
    x <- matrix(x, ncol = 1)
    if (!inherits(x, "dgCMatrix")) x <- as(x, "dgCMatrix")
    return(x)
  })
  max_rows <- max(vapply(mats, nrow, numeric(1)))
  mats_padded <- lapply(mats, function(x) {
    r <- nrow(x)
    if (r < max_rows) {
      padding <- Matrix(0, nrow = max_rows - r, ncol = ncol(x), sparse = TRUE)
      x <- rbind2(x, padding)
    }
    return(x)
  })
  result <- do.call(cbind, mats_padded)
  return(result)
}

calc_rmse <- function(x, y) sqrt(mean((x - y)^2))

index_spliter <- function(array, n_folds = 5){
  len <- length(array)
  fold_id = sample(rep(1:n_folds, length.out = len))
  split(seq_len(len), fold_id)
}

# sieve -------------------------------------------------------------------
# Motivated from the codes by Zhang, T., & Simon, N. (2023).
# Regression in tensor product spaces by the method of sieves. 
# Electronic journal of statistics, 17(2), 3660.


#' Cross-Validated High-Dimensional Sieve Estimation
#'
#' @description
#' Performs $K$-fold cross-validation to select the optimal number of sieve bases ($J$) 
#' and the $L_1$ penalty parameter ($\lambda$) for high-dimensional sieve regression.
#' This function is optimized for computational efficiency and strict statistical validity, 
#' specifically designed for tensor product bases (e.g., Cosine basis) with nested structures.
#'
#' @details
#' The algorithm incorporates two major optimizations to handle high-dimensional settings efficiently:
#' 
#' 1. **Global Penalty Sequence Pre-calculation**: To strictly avoid data leakage and maintain 
#'    the symmetry (exchangeability) of cross-validation, the global $\lambda$ sequence is 
#'    pre-calculated using the full dataset for each candidate basis size. This guarantees that 
#'    the validation errors across different folds are evaluated on the exact same penalty grid.
#' 
#' 2. **Loop Inversion & Nested Subspace Slicing**: To minimize the computational bottleneck of 
#'    evaluating complex basis functions (e.g., trigonometric transformations), the function 
#'    applies `sieve_preprocess` only once per CV fold using the maximum basis size ($J_{max}$). 
#'    For smaller candidate sizes, it utilizes the nested property of the basis space by 
#'    directly slicing the pre-computed design matrix (`Phi`) and `index_matrix` in $O(1)$ time.
#'
#' @param X A numeric matrix of predictors. The first column is assumed to be an intercept or ID.
#' @param y A numeric vector of the response variable.
#' @param basis_numbers A numeric vector specifying the candidate numbers of sieve bases to tune. 
#'        If NULL, a default heuristic sequence based on the sample size and dimension is used.
#' @param n_folds An integer specifying the number of cross-validation folds (default: 5).
#' @param type A character string specifying the type of sieve basis (default: 'cosine').
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{best_lambda}: The optimal $L_1$ penalty value minimizing the CV MSE.
#'   \item \code{best_basis_number}: The optimal number of sieve bases ($J$).
#'   \item \code{cv_results}: A \code{data.table} summarizing the average validation MSE across the parameter grid.
#' }
cv.sieve <- function(X, y, basis_numbers = NULL, n_folds = 5, type = 'cosine') {
  n <- nrow(X)
  xdim <- ncol(X)
  if (is.null(basis_numbers)) {
    basis_numbers = ceiling(
      c(xdim * c(5, n^(1 / 5), n^(1 / 3), n^(1 / 2))
        , xdim^2
        # , xdim^2 * c(5, n^(1 / 5), n^(1 / 3))
        )
    )
  }
  basis_numbers <- sort(unique(basis_numbers)) 
  max_basisN <- max(basis_numbers)
  validation_split_index <- index_spliter(1:n, n_folds = n_folds)
  # ------------------------------------------------------------------------
  # [Phase 1] Global Lambda Sequence Pre-calculation
  # ------------------------------------------------------------------------
  global_lambda_list <- list()
  full_model <- sieve_preprocess(X, max_basisN, type = type)
  for (basisN in basis_numbers) {
    # Avoid deep copy: only update necessary fields
    current_model <- full_model
    current_model$Phi <- full_model$Phi[, 1:basisN, drop = FALSE]
    current_model$index_matrix <- full_model$index_matrix[1:basisN, , drop = FALSE]
    current_model$basisN <- basisN
    dummy_fit <- sieve_solver(current_model, y)
    global_lambda_list[[as.character(basisN)]] <- dummy_fit$lambda
  }
  # ------------------------------------------------------------------------
  # [Phase 2] Cross-Validation (preallocate results)
  # ------------------------------------------------------------------------
  result_list <- vector("list", n_folds * length(basis_numbers))
  res_idx <- 1
  for (split_index in seq_len(n_folds)) {
    train_idx <- setdiff(1:n, validation_split_index[[split_index]])
    valid_idx <- validation_split_index[[split_index]]
    Train_X <- X[train_idx, , drop = FALSE]
    Train_Y <- y[train_idx]
    Valid_X <- X[valid_idx, , drop = FALSE]
    Valid_Y <- y[valid_idx]
    full_train_model <- sieve_preprocess(Train_X, max_basisN, type = type)
    for (basisN in basis_numbers) {
      current_train_model <- full_train_model
      current_train_model$Phi <- full_train_model$Phi[, 1:basisN, drop = FALSE]
      current_train_model$index_matrix <- full_train_model$index_matrix[1:basisN, , drop = FALSE]
      current_train_model$basisN <- basisN
      current_lambda_seq <- global_lambda_list[[as.character(basisN)]]
      sieve.fit <- sieve_solver(current_train_model, Train_Y, lambda = current_lambda_seq)
      sieve.validation <- sieve_predict(model = sieve.fit, testX = Valid_X, testY = Valid_Y)
      result_list[[res_idx]] <- data.table(
        l1_penalty = current_lambda_seq,
        l1_penalty_index = seq_along(current_lambda_seq),
        basisN = rep(basisN, length(current_lambda_seq)),
        validation_MSE = sieve.validation$MSE,
        split_index = rep(split_index, length(current_lambda_seq))
      )
      res_idx <- res_idx + 1
    }
  }
  parameter_tuning_reference <- rbindlist(result_list)
  # ------------------------------------------------------------------------
  # [Phase 3] Aggregation & Best Parameter Selection
  # ------------------------------------------------------------------------
  average_data <- parameter_tuning_reference[, .(
    average_MSE = mean(validation_MSE),
    l1_penalty = mean(l1_penalty) 
  ), by = .(l1_penalty_index, basisN)]
  best_combination_index <- which.min(average_data$average_MSE)
  best_lambda <- average_data$l1_penalty[best_combination_index]
  best_basis_number <- average_data$basisN[best_combination_index]
  return(list(
    best_lambda = best_lambda,
    best_basis_number = best_basis_number,
    cv_results = average_data 
  ))
}

cv_fitH_from_prxy <- function(X, Z, h = "sieve", n_folds = 5) {
  p1 <- ncol(X)
  p2 <- ncol(Z)
  n <- nrow(X)
  n_cores <- parallel::detectCores() - 2
  if (h == "linear") {
    # Linear fmap with intercept (no centering assumed)
    fit.list <- pbmcapply::pbmclapply(1:p2, function(j) {
      cv.glmnet(X, Z[, j], alpha = 1)
    }, mc.cores = n_cores)
    ϴ <- lapply(fit.list, function(fit)
      coef(fit, s = 'lambda.min'))  # keep intercept row
    return(do.call(cbind, ϴ))
  }
  # Sieve
  sv <- pbmcapply::pbmclapply(1:p2, function(j) 
    cv.sieve(X, Z[, j]), mc.cores = n_cores)
  sv <- pbmcapply::pbmclapply(1:p2, function(j) {
    sv_prep <- sieve_preprocess(X, sv[[j]]$best_basis_number, type = 'cosine')
    sv_fit <- sieve_solver(sv_prep, Z[, j], lambda = sv[[j]]$best_lambda)
    as.numeric(sv_fit$beta_hat)  # keep full beta_hat: Phi[,1]=1 carries the intercept
  }, mc.cores = n_cores)
  ϴ <- sv %>% do.call(padder0, .)
  ϴ
}