htl.train <- function(p_coefs, X, Y, Z = NULL) {
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
      proxy_fit <- cv.glmnet(D, proxyk$y, alpha = 1)
      list(p_coefs = coef(proxy_fit, s = 'lambda.1se'), h_Θ = matrix(0, 2, 2))
    } else {
      D <- cbind(proxyk$D$X, proxyk$D$Z)
      proxy_fit <- cv.glmnet(D, proxyk$y, alpha = 1)
      h_Θ <- cv.fmap(proxyk$D$X, proxyk$D$Z, h = fmap)
      list(p_coefs = coef(proxy_fit, s = 'lambda.1se'), h_Θ = h_Θ)
    }
  })
  proxy_h_Θs <- lapply(proxy_fits, `[[`, "h_Θ")
  proxy_h_Θ_avg <- sum_padded_sparse_matrices(proxy_h_Θs) / K
  proxy_coefs_mat <- lapply(proxy_fits, `[[`, "p_coefs") %>% do.call(cbind, .)
  proxy_coef_avg <- rowMeans(as.matrix(proxy_coefs_mat)) 
  list(p_coefs = proxy_coef_avg, h_Θ = proxy_h_Θ_avg, fmap_type = fmap)
}

htl <- function(proxy_fit, train, test) {
  fmap_type <- proxy_fit$fmap_type
  if (fmap_type == "none") {
    # no imputation
    test_D <- test$X
    p_coefs <- proxy_fit$p_coefs[1:(ncol(test_D) + 1)]
    z_trgt <- NULL
  } else {
    p_coefs <- proxy_fit$p_coefs
    z_trn <- cv.fmap.predict(proxy_fit$h_Θ, train$X, fmap_type)
    z_test <- cv.fmap.predict(proxy_fit$h_Θ, test$X, fmap_type)
    test_D <- cbind(test$X, z_test)
  }
  # Debiasing on train data
  train_tl_fit <- htl.train(p_coefs, train$X, train$y, z_trn)
  # Coefficient combination
  h_d_t <- coef(train_tl_fit, s = 'lambda.1se')
  trn_coef <- as.vector(p_coefs) + as.vector(h_d_t)
  h_y_test <- cbind(1, test_D) %*% trn_coef
  pe <- .calc_rmse(h_y_test, test$y)
  return(list(pe = pe, coef = trn_coef))
}

sim.runner <- function(seedno = 1, K = 1, np = 500, nt = 30, ntest = 100,
                       p1 = 30, p2 = 30, dgp_fmap = 'nlinear', sp = 'nsp') {
  cnfg1 <- dgp.confg(seedno, K, p1, p2, dgp_fmap)
  D_data <- dgp.gen.D(seedno, dgp_fmap, np, nt, ntest, cnfg1$Prxy, cnfg1$Trgt)
  
  train <- D_data$traindata
  test   <- D_data$testdata
  beta_t <- D_data$βt
  proxy_list <- D_data$proxydata[[sp]]
  
  res_lasso  <- vanilla.glmnet(train, test, beta_t, 1)
  
  proxy_fit_nl <- agg.proxy(proxy_list, 'nlinear')
  res_htl_nl <- htl(proxy_fit_nl, train, test, beta_t)
  
  proxy_fit_ln <- agg.proxy(proxy_list, 'linear')
  res_htl_ln <- htl(proxy_fit_ln, train, test, beta_t)
  
  proxy_fit_none <- agg.proxy(proxy_list, 'none')
  res_hmtl <- htl(proxy_fit_none, train, test, beta_t)
  
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