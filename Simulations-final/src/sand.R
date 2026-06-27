# =============================================================================
# SAND: Signal-Based Adaptive Negative Transfer Defense
# Linear feature mapping case
# =============================================================================
make_proxy_preds_linear <- function(target, proxy_list) {
  K  <- length(proxy_list)
  X_t <- as.matrix(target$X)          # n_t × p1
  nt <- nrow(X_t)

  lapply(seq_len(K), function(k) {
    prxk <- proxy_list[[k]]
    Xp   <- as.matrix(prxk$D$X)      # n_p × p1
    Zp   <- as.matrix(prxk$D$Z)      # n_p × p2
    Yp   <- as.numeric(prxk$Y)       # n_p

    # (a) Linear feature map estimation: hat(P)_p via column-wise Lasso
    P_hat <- cv.fmap(Xp, Zp, 'linear')

    # (b) Proxy regression: fit OLS on full proxy design [X_p, Z_p]
    Dp    <- cbind(Xp, Zp)           # n_p × p
    lm_fit <- lm(Yp ~ Dp)
    omega_hat <- coef(lm_fit)        # length p+1 (intercept + p coefs)

    # Target marginalized prediction:
    p1        <- ncol(Xp)
    p2        <- ncol(Zp)
    omega1    <- omega_hat[1:(p1 + 1)]
    omega2    <- omega_hat[(p1 + 2):(p1 + p2 + 1)]

    # Predicted target response
    as.numeric(cbind(1, X_t) %*% (omega1 + P_hat %*% omega2))
  })
}


#' Signal-Based Adaptive Negative Transfer Defense (SAND)
run_sand <- function(X_t, y_t, proxy_preds_list, gamma_scale = 1) {
  X_t <- as.matrix(X_t)
  y_t <- as.numeric(y_t)
  nt  <- nrow(X_t)
  p1  <- ncol(X_t)
  K   <- length(proxy_preds_list)

  # ==================================================================
  # Step 1: Compute Target Vanilla Baseline
  # ==================================================================
  cv_base  <- cv.glmnet(X_t, y_t, alpha = 1)
  idx_1se <- which.min(abs(cv_base$lambda - cv_base$lambda.1se))
  S_base  <- cv_base$cvm[idx_1se]
  
  idx_min <- which(cv_base$lambda == cv_base$lambda.min)
  s_hat   <- cv_base$nzero[idx_min]
  y_hat   <- predict(cv_base, newx = X_t, s = "lambda.min")
  rss     <- sum((y_t - drop(y_hat))^2)
  df_adj     <- max(nt - s_hat, 1) 
  tau_sq_hat <- rss / df_adj

  gamma_n <- gamma_scale * (tau_sq_hat + S_base) * sqrt(log(p1) / nt)
  # ==================================================================
  # Step 2 & 3: Evaluate HTL Proxy Residuals & Adaptive Thresholding
  # ==================================================================
  tS_vec <- sapply(seq_len(K), function(k) {
    y_hat_k <- as.numeric(proxy_preds_list[[k]])
    mean((y_t - y_hat_k)^2)
  })

  # Hard Thresholding
  threshold <- S_base + gamma_n
  hA        <- which(tS_vec <= threshold)

  if (length(hA) == 0L) {
    return(integer(0))
  }
  return(as.integer(hA))
}


#' Full SAND pipeline: DGP → Proxy Predictions → SAND selection
run_sand_pipeline <- function(D_data, sp = "nsp", gamma_n = NULL) {
  target      <- D_data$targetdata
  proxy_list  <- D_data$proxydata[[sp]]

  proxy_preds <- make_proxy_preds_linear(target, proxy_list)

  hA <- run_sand(
    X_t              = target$X,
    y_t              = target$Y,
    proxy_preds_list = proxy_preds,
    gamma_n          = gamma_n
  )
  hA
}