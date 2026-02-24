# Simulation: Heterogeneous Transfer Learning vs Homogeneous Transfer Learning vs Lasso
# Compares HTL (with linear and sieve feature imputation) vs HmTL vs Lasso-only baseline

library(glmnet)
library(tidyverse)
library(data.table)

source('../sv-cv.R')
source('../sv-prim.R')

# ==============================================================================
# Helper functions
# ==============================================================================

calc_rmse <- function(y_true, y_pred) {
  sqrt(mean((y_true - y_pred)^2))
}

calc_mae <- function(y_true, y_pred) {
  mean(abs(y_true - y_pred))
}

# ==============================================================================
# Data Generation Functions
# ==============================================================================

# Generate source data with more features than target
generate_source_data <- function(n, p1, regime = 'additive') {
  X <- matrix(rnorm(n * p1), nrow = n, ncol = p1)
  
  if (regime == 'additive') {
    # Simple additive model
    Y <- rowSums(X[, 1:5]) + rnorm(n, sd = 0.5)
  } else if (regime == 'multiplicative') {
    # Multiplicative interactions
    Y <- X[, 1] * X[, 2] + X[, 3] + rnorm(n, sd = 1)
  } else if (regime == 'nonlinear') {
    # Non-linear transformations
    Y <- sin(X[, 1]) + cos(X[, 2]) + X[, 3]^2 + rnorm(n, sd = 0.5)
  }
  
  colnames(X) <- paste0('X_s', 1:p1)
  list(X = X, Y = Y, p1 = p1)
}

# Generate target data with subset of source features
generate_target_data <- function(n, p2, p1_full, regime = 'additive') {
  # Select p2 features from the p1 available in source
  selected_features <- sample(1:p1_full, size = p2, replace = FALSE)
  X_full <- matrix(rnorm(n * p1_full), nrow = n, ncol = p1_full)
  X <- X_full[, selected_features, drop = FALSE]
  
  if (regime == 'additive') {
    Y <- rowSums(X_full[, 1:min(5, p1_full)]) + rnorm(n, sd = 0.5)
  } else if (regime == 'multiplicative') {
    Y <- X_full[, 1] * X_full[, 2] + X_full[, 3] + rnorm(n, sd = 1)
  } else if (regime == 'nonlinear') {
    Y <- sin(X_full[, 1]) + cos(X_full[, 2]) + X_full[, 3]^2 + rnorm(n, sd = 0.5)
  }
  
  colnames(X) <- paste0('X_t', 1:p2)
  list(X = X, Y = Y, X_full = X_full, selected_features = selected_features, p2 = p2)
}

# Generate homogeneous data (same features in source and target)
generate_homogeneous_data <- function(n, p, regime = 'additive') {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  if (regime == 'additive') {
    Y <- rowSums(X[, 1:min(5, p)]) + rnorm(n, sd = 0.5)
  } else if (regime == 'multiplicative') {
    Y <- X[, 1] * X[, 2] + X[, 3] + rnorm(n, sd = 1)
  } else if (regime == 'nonlinear') {
    Y <- sin(X[, 1]) + cos(X[, 2]) + X[, 3]^2 + rnorm(n, sd = 0.5)
  }
  
  colnames(X) <- paste0('X', 1:p)
  list(X = X, Y = Y, p = p)
}

# ==============================================================================
# Method 1: HTL with Linear Mapping (maps source common features to target features)
# ==============================================================================

htl_linear <- function(X_source, Y_source, X_target_train, Y_target_train, 
                       X_target_test, Y_target_test) {
  
  # Assume the target features are the first p2 features of source
  p2 <- ncol(X_target_train)
  p1 <- ncol(X_source)
  
  # Fit linear mappings from source common features to target features
  # using the first p2 features of source data
  X_source_common <- X_source[, 1:p2, drop = FALSE]
  X_source_to_map <- X_source[, (p2+1):p1, drop = FALSE]
  
  # If there are features to map
  if (ncol(X_source_to_map) > 0) {
    fmap <- cv_fitH_from_prxy(X_source_common, X_source_to_map, 
                               type = "linear", n_folds = 5)
    
    # Predict mapped features on target training set
    X_mapped_train <- fmap_predict(fmap, X_target_train, type = 'linear')
    X_combined_train <- cbind(X_target_train, X_mapped_train)
    
    # Predict mapped features on target test set
    X_mapped_test <- fmap_predict(fmap, X_target_test, type = 'linear')
    X_combined_test <- cbind(X_target_test, X_mapped_test)
  } else {
    X_combined_train <- X_target_train
    X_combined_test <- X_target_test
  }
  
  # Fit lasso on combined features (target + mapped proxy features)
  lasso_fit <- cv.glmnet(X_combined_train, Y_target_train, alpha = 1, nfolds = 5)
  
  # Predict on test set
  Y_pred <- predict(lasso_fit, X_combined_test, s = 'lambda.1se')
  
  list(
    method = "HTL_Linear",
    Y_pred = as.numeric(Y_pred),
    rmse = calc_rmse(Y_target_test, Y_pred),
    mae = calc_mae(Y_target_test, Y_pred),
    lasso_fit = lasso_fit,
    fmap = if (ncol(X_source_to_map) > 0) fmap else NULL
  )
}

# ==============================================================================
# Method 2: HTL with Sieve Estimation (non-linear mapping via sieves)
# ==============================================================================

htl_sieve <- function(X_source, Y_source, X_target_train, Y_target_train, 
                      X_target_test, Y_target_test) {
  
  p2 <- ncol(X_target_train)
  p1 <- ncol(X_source)
  
  # Fit sieve mappings from source common features to target features
  X_source_common <- X_source[, 1:p2, drop = FALSE]
  X_source_to_map <- X_source[, (p2+1):p1, drop = FALSE]
  
  # If there are features to map
  if (ncol(X_source_to_map) > 0) {
    fmap <- cv_fitH_from_prxy(X_source_common, X_source_to_map, 
                               type = "sieve", n_folds = 5)
    
    # Predict mapped features on target training set using sieve
    X_mapped_train <- fmap_predict(fmap, X_target_train, type = 'sieve')
    X_combined_train <- cbind(X_target_train, X_mapped_train)
    
    # Predict mapped features on target test set using sieve
    X_mapped_test <- fmap_predict(fmap, X_target_test, type = 'sieve')
    X_combined_test <- cbind(X_target_test, X_mapped_test)
  } else {
    X_combined_train <- X_target_train
    X_combined_test <- X_target_test
  }
  
  # Fit lasso on combined features (target + mapped proxy features)
  lasso_fit <- cv.glmnet(X_combined_train, Y_target_train, alpha = 1, nfolds = 5)
  
  # Predict on test set
  Y_pred <- predict(lasso_fit, X_combined_test, s = 'lambda.1se')
  
  list(
    method = "HTL_Sieve",
    Y_pred = as.numeric(Y_pred),
    rmse = calc_rmse(Y_target_test, Y_pred),
    mae = calc_mae(Y_target_test, Y_pred),
    lasso_fit = lasso_fit,
    fmap = if (ncol(X_source_to_map) > 0) fmap else NULL
  )
}

# ==============================================================================
# Method 3: Homogeneous Transfer Learning (HmTL)
# Uses source data with same features as target
# ==============================================================================

hmtl <- function(X_source, Y_source, X_target_train, Y_target_train, 
                 X_target_test, Y_target_test) {
  
  # Combine source and target training data
  X_combined <- rbind(X_source, X_target_train)
  Y_combined <- c(Y_source, Y_target_train)
  
  # Fit lasso on combined data
  lasso_fit <- cv.glmnet(X_combined, Y_combined, alpha = 1, nfolds = 5)
  
  # Predict on target test set
  Y_pred <- predict(lasso_fit, X_target_test, s = 'lambda.1se')
  
  list(
    method = "HmTL",
    Y_pred = as.numeric(Y_pred),
    rmse = calc_rmse(Y_target_test, Y_pred),
    mae = calc_mae(Y_target_test, Y_pred),
    lasso_fit = lasso_fit
  )
}

# ==============================================================================
# Method 4: Lasso on Target Data Only (baseline)
# ==============================================================================

lasso_target_only <- function(X_target_train, Y_target_train, 
                              X_target_test, Y_target_test) {
  
  # Fit lasso only on target training data
  lasso_fit <- cv.glmnet(X_target_train, Y_target_train, alpha = 1, nfolds = 5)
  
  # Predict on target test set
  Y_pred <- predict(lasso_fit, X_target_test, s = 'lambda.1se')
  
  list(
    method = "Lasso_Target_Only",
    Y_pred = as.numeric(Y_pred),
    rmse = calc_rmse(Y_target_test, Y_pred),
    mae = calc_mae(Y_target_test, Y_pred),
    lasso_fit = lasso_fit
  )
}

# ==============================================================================
# Main Simulation
# ==============================================================================

main <- function(n_sim = 100, 
                 n_source = 500, 
                 n_target_train = 100, 
                 n_target_test = 100,
                 p1 = 15,      # number of source features
                 p2 = 8,       # number of target features
                 regime = 'additive') {
  
  set.seed(42)
  results_list <- list()
  
  cat(sprintf("Starting simulation with %d iterations\n", n_sim))
  cat(sprintf("Source: %d samples, %d features\n", n_source, p1))
  cat(sprintf("Target: %d train + %d test samples, %d features\n", 
              n_target_train, n_target_test, p2))
  cat(sprintf("Regime: %s\n\n", regime))
  
  for (sim_id in 1:n_sim) {
    if (sim_id %% 10 == 0) {
      cat(sprintf("Completed %d/%d simulations\n", sim_id, n_sim))
    }
    
    # Generate source data (heterogeneous - has extra features)
    source_data <- generate_source_data(n_source, p1, regime)
    
    # Generate target data (heterogeneous - subset of source features)
    target_data <- generate_target_data(n_target_train + n_target_test, 
                                        p2, p1, regime)
    X_target_all <- target_data$X
    Y_target_all <- target_data$Y
    
    X_target_train <- X_target_all[1:n_target_train, , drop = FALSE]
    Y_target_train <- Y_target_all[1:n_target_train]
    X_target_test <- X_target_all[(n_target_train + 1):(n_target_train + n_target_test), , drop = FALSE]
    Y_target_test <- Y_target_all[(n_target_train + 1):(n_target_train + n_target_test)]
    
    # Generate homogeneous source data for HmTL
    source_homo <- generate_homogeneous_data(n_source, p2, regime)
    
    # Run methods
    tryCatch({
      result_htl_linear <- htl_linear(source_data$X, source_data$Y,
                                      X_target_train, Y_target_train,
                                      X_target_test, Y_target_test)
      results_list[[paste0(sim_id, "_HTL_Linear")]] <- result_htl_linear
    }, error = function(e) {
      cat("Error in HTL_Linear:", conditionMessage(e), "\n")
    })
    
    tryCatch({
      result_htl_sieve <- htl_sieve(source_data$X, source_data$Y,
                                    X_target_train, Y_target_train,
                                    X_target_test, Y_target_test)
      results_list[[paste0(sim_id, "_HTL_Sieve")]] <- result_htl_sieve
    }, error = function(e) {
      cat("Error in HTL_Sieve:", conditionMessage(e), "\n")
    })
    
    tryCatch({
      result_hmtl <- hmtl(source_homo$X, source_homo$Y,
                          X_target_train, Y_target_train,
                          X_target_test, Y_target_test)
      results_list[[paste0(sim_id, "_HmTL")]] <- result_hmtl
    }, error = function(e) {
      cat("Error in HmTL:", conditionMessage(e), "\n")
    })
    
    tryCatch({
      result_lasso <- lasso_target_only(X_target_train, Y_target_train,
                                        X_target_test, Y_target_test)
      results_list[[paste0(sim_id, "_Lasso")]] <- result_lasso
    }, error = function(e) {
      cat("Error in Lasso:", conditionMessage(e), "\n")
    })
  }
  
  cat("\nSimulation complete. Summarizing results...\n")
  
  # Summarize results
  summary_dt <- do.call(rbind, lapply(results_list, function(res) {
    data.table(
      method = res$method,
      rmse = res$rmse,
      mae = res$mae
    )
  }))
  
  # Calculate summary statistics by method
  summary_stats <- summary_dt[, .(
    mean_rmse = mean(rmse),
    sd_rmse = sd(rmse),
    median_rmse = median(rmse),
    mean_mae = mean(mae),
    sd_mae = sd(mae),
    median_mae = median(mae),
    n = .N
  ), by = method]
  
  list(
    results = results_list,
    summary = summary_dt,
    summary_stats = summary_stats
  )
}

# ==============================================================================
# Execution (uncomment to run)
# ==============================================================================

# results <- main(n_sim = 50, n_source = 500, n_target_train = 100, 
#                 n_target_test = 100, p1 = 15, p2 = 8, regime = 'additive')
# print(results$summary_stats)
