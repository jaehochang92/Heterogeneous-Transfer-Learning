htl.target <- function(p_coefs, X, Y, Z = NULL) {
  if (is.null(Z)) {
    Dt <- X
  } else {
    Dt <- cbind(X, Z)
  }
  Dt <- as.matrix(Dt)
  p_coefs <- p_coefs[1:(ncol(Dt) + 1)]
  cv.glmnet(Dt, Y, alpha = 1, offset = cbind(1, Dt) %*% p_coefs)
}

agg.proxy <- function(proxy_list, fmap = c('none', 'linear', 'nlinear')) {
  fmap <- match.arg(fmap)
  K <- length(proxy_list)
  proxy_fits <- lapply(proxy_list, function(proxyk) {
    if (fmap == 'none') {
      D <- proxyk$D$X
      proxy_fit <- lm(proxyk$Y ~ D)
      list(p_coefs = coef(proxy_fit), h_Θ = matrix(0, 2, 2))
    } else {
      D <- cbind(proxyk$D$X, proxyk$D$Z)
      proxy_fit <- lm(proxyk$Y ~ D)
      h_Θ <- cv.fmap(proxyk$D$X, proxyk$D$Z, h = fmap)
      list(p_coefs = coef(proxy_fit), h_Θ = h_Θ)
    }
  })
  proxy_h_Θs <- lapply(proxy_fits, `[[`, "h_Θ")
  proxy_h_Θ_avg <- sum_padded_sparse_matrices(proxy_h_Θs) / K
  proxy_coefs_mat <- lapply(proxy_fits, `[[`, "p_coefs") %>% do.call(cbind, .)
  proxy_coef_avg <- rowMeans(as.matrix(proxy_coefs_mat)) 
  list(p_coefs = proxy_coef_avg, h_Θ = proxy_h_Θ_avg, fmap_type = fmap)
}

htl.eval <- function(proxy_fit, target, test, beta_true) {
  fmap_type <- proxy_fit$fmap_type
  if (fmap_type == "none") {
    # no imputation
    test_D <- test$X
    p_coefs <- proxy_fit$p_coefs[1:(ncol(test_D) + 1)]
    z_trgt <- NULL
    eval_beta <- beta_true[1:ncol(target$X)]
  } else {
    p_coefs <- proxy_fit$p_coefs
    z_trgt <- cv.fmap.predict(proxy_fit$h_Θ, target$X, fmap_type)
    z_test <- cv.fmap.predict(proxy_fit$h_Θ, test$X, fmap_type)
    test_D <- cbind(test$X, z_test)
    eval_beta <- beta_true
  }
  # Debiasing on target data
  target_tl_fit <- htl.target(p_coefs, target$X, target$Y, z_trgt)
  # Coefficient combination
  h_d_t <- coef(target_tl_fit, s = 'lambda.min')
  trgt_coef <- as.vector(p_coefs) + as.vector(h_d_t)
  est_err <- .calc_rmse(trgt_coef[-1], eval_beta)
  h_y_test <- cbind(1, test_D) %*% trgt_coef
  pe <- .calc_rmse(h_y_test, test$Y)
  return(list(est_err = est_err, pe = pe, coef = trgt_coef))
}

vanilla.glmnet <- function(target, test, beta_true, alp) {
  lasso_fit <- cv.glmnet(as.matrix(target$X), target$Y, alpha = alp)
  lasso_coef <- as.vector(coef(lasso_fit, s = 'lambda.1se'))
  est_err <- .calc_rmse(lasso_coef[-1], beta_true[1:(length(lasso_coef) - 1)])
  h_y_test <- predict(lasso_fit, test$X, s = 'lambda.1se')
  pe <- .calc_rmse(h_y_test, test$Y)
  return(list(est_err = est_err, pe = pe))
}

sim.runner <- function(seedno = 1, K = 1, np = 500, nt = 30, ntest = 100,
                       p1 = 30, p2 = 30, dgp_fmap = 'nlinear', sp = 'nsp') {
  cnfg1 <- dgp.confg(seedno, K, p1, p2, dgp_fmap)
  D_data <- dgp.gen.D(seedno, dgp_fmap, np, nt, ntest, cnfg1$Prxy, cnfg1$Trgt)
  
  target <- D_data$targetdata
  test   <- D_data$testdata
  beta_t <- D_data$βt
  proxy_list <- D_data$proxydata[[sp]]
  
  res_lasso  <- vanilla.glmnet(target, test, beta_t, 1)
  
  proxy_fit_nl <- agg.proxy(proxy_list, 'nlinear')
  res_htl_nl <- htl.eval(proxy_fit_nl, target, test, beta_t)
  
  proxy_fit_ln <- agg.proxy(proxy_list, 'linear')
  res_htl_ln <- htl.eval(proxy_fit_ln, target, test, beta_t)
  
  proxy_fit_none <- agg.proxy(proxy_list, 'none')
  res_hmtl <- htl.eval(proxy_fit_none, target, test, beta_t)
  
  est_err_vec <- c(Lasso = res_lasso$est_err, 
                   HmTL = res_hmtl$est_err, 
                   HTL.LN = res_htl_ln$est_err, 
                   HTL.NL = res_htl_nl$est_err)
  
  pe_vec <- c(Lasso = res_lasso$pe, 
              HmTL = res_hmtl$pe, 
              HTL.LN = res_htl_ln$pe, 
              HTL.NL = res_htl_nl$pe,
              Oracle = sqrt(mean(D_data$testdata$εtest^2)))
  
  return(list(Est_Error = est_err_vec, Prediction_Error = pe_vec))
}