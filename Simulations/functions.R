# Transfer Learning with missing covariates - sieve
source('sv-prim.R')

library(MASS)
library(glmnet)
library(stargazer)
library(pbapply)
library(tidyr)
library(truncnorm)
library(Sieve)

# https://github.com/terrytianyuzhang/Sieve/blob/master/R/SieveFittingModels.R

# DGP ---------------------------------------------------------------------


GenD <- function(n, p1, p2, Σx, Σξ, A, alg = 'jae', trgt = F) {
  if (alg == 'jae') {
    if (!is.null(Σx)) {
      X <- mvrnorm(n, rep(0, p1), Σx)
    } else {
      X <- matrix(runif(n * p1, -2, 2), n)
    }
    tb <- tibble(X = X)
    if (p2 > 0) {
      Ξ <- mvrnorm(n, rep(0, p2), Σξ)
      Z <- H(X, p2, A, trgt) + Ξ
      tb$Z <- Z
      # pairs(data.frame(X = X, Z = Z), cex = .1, main = n)
      # Sys.sleep(1)
    }
  }
  return(tb)
}

### Data-Generating Process
# K: # of proxy studies
# np: sample size of proxy studies (> d, set to be identical across K studies)
# nt: sample size of target study (typically nt{\gg}log(p1 + p2))
# p1, p2: # of matched and mismatched features (may depend on nt)
# sparse: whether beta contrast should be sparse (proximity)
DGP_transfer = function(K, np, nt, p1, p2, ntest, Prxy, Trgt) {
  p <- p1 + p2
  𝒫 <- list()
  δst <- list()
  E.cum <- mvrnorm(max(np * K, nt + ntest), c(0, 0), matrix(c(1, .5, .5, 1), 2))
  A <- replicate(p2, 1:5)
  # A <- replicate(p2, 1:round(log(p1, 2)))
  # A <- replicate(p2, sample(1:p1, size = round(log(p1, 2))))
  
  ## Proxy domain
  for (k in 1:K) {
    εp <- E.cum[(np * (k - 1) + 1):(np * k), 2]
    ## difference of beta1
    Dp <- GenD(np, p1, p2, Prxy$Σxp[[k]][1:p1, 1:p1], 
               Prxy$Σξp[[k]][1:p2, 1:p2], A)
    δst[['sp']] <- Prxy$δst_sprs[[k]][1:p]
    δst[['nsp']] <- Prxy$δst_nsprs[[k]][1:p]
    for (δsti in c('sp', 'nsp')) {
      βpk <- Trgt$βt[1:p] - δst[[δsti]][1:p]
      Yp <- as.matrix(Dp) %*% βpk + εp
      proxydata = tibble(Y = c(Yp), D = Dp)
      𝒫[[δsti]][[k]] <- proxydata
    }
  }
  
  ## Target domain
  εt     = E.cum[1:nt, 1]
  εtest  = E.cum[1:ntest + nt, 1]
  Dtrgt <- GenD(nt + ntest, p1, p2, NULL, Trgt$Σξt[1:p2, 1:p2], A, trgt = T)
  Dt <- Dtrgt[1:nt,]
  Dtest <- Dtrgt[1:ntest + nt,]
  Yt     = as.matrix(Dt) %*% Trgt$βt[1:p] + εt
  Ytest  = as.matrix(Dtest) %*% Trgt$βt[1:p] + εtest
  
  𝒯 <- tibble(Y = c(Yt), X = Dt$X) # missing Z
  testdata = tibble(Y = c(Ytest), D = Dtest) # keep Z for oracle error
  
  return(
    list(
      proxydata = 𝒫,
      targetdata = 𝒯,
      testdata = testdata,
      p1 = p1,
      p2 = p2,
      K = K,
      βt = Trgt$βt[1:p],
      Θt = Trgt$Θt[1:p1, 1:p2]
    )
  )
}


# Estimation --------------------------------------------------------------


### Method 1: Target Lasso only

lasso_trgt = function(𝒟) {
  p1 = 𝒟$p1
  Yt = 𝒟$targetdata$Y
  Xt = 𝒟$targetdata$X
  Ytest = 𝒟$testdata$Y
  Xtest = 𝒟$testdata$D$X
  
  lassotarget <- cv.glmnet(Xt, Yt, alpha = 1, intercept = F)
  hβ1t <- coef(lassotarget, 'lambda.min')
  Ypred <- predict(lassotarget, newx = Xtest)
  pe = norm(Ytest - Ypred, '2')
  ee = norm(𝒟$βt[1:p1] - hβ1t[1:p1 + 1], '2')
  return(c(pe = pe, ee = ee))
}

### Method 2: Homogeneous TL (Bastani)

TL_matched = function(𝒟, δsti) {
  p1 = 𝒟$p1
  hβp1 <- rep(0, p1)
  for (k in 1:𝒟$K) {
    Prxyk <- 𝒟$proxydata[[δsti]][[k]]
    Yp  = Prxyk$Y
    Xp = Prxyk$D$X
    estproxy = lm(Yp ~ Xp - 1)
    hβp1 <- hβp1 + coef(estproxy) / 𝒟$K
  }
  Xt = 𝒟$targetdata$X
  Xtest = 𝒟$testdata$D$X
  Yt = 𝒟$targetdata$Y
  Ytest = 𝒟$testdata$Y
  
  tYt         = c(Yt - Xt %*% hβp1)
  combinedlasso = cv.glmnet(Xt, tYt, alpha = 1, intercept = F)
  hdelta        = coef(combinedlasso)
  Ypred         = predict(combinedlasso, newx = Xtest) + c(Xtest %*% hβp1)
  pe = norm(Ytest - Ypred, '2')
  ee = norm(𝒟$βt[1:p1] - hdelta[1:p1 + 1] - hβp1, '2')
  return(c(pe = pe, ee = ee))
}


### Method 4 (Proposed): Two-Stage Joint Estimator - including the imputation stage which uses X as a predictor


TL_impute = function(𝒟, δsti) {
  p1 = 𝒟$p1
  p2 = 𝒟$p2
  p <- p1 + p2
  hβp <- rep(0, p)
  M <- 10 * p
  hΘ <- matrix(0)
  for (k in 1:𝒟$K) {
    hΘk <- 0
    Prxyk <- 𝒟$proxydata[[δsti]][[k]]
    Yp  = Prxyk$Y
    Xp = Prxyk$D$X
    Zp = Prxyk$D$Z
    for (j in 1:p2) {
      df <- data.frame(Zp[, j], Xp)
      colnames(df) <- c('Y', paste0('X', 1:p1))
      sv.cv.fit <- jae_sv_cv(df, type = 'cosine', n_folds = 5)
      hΘk <- cbind_pad0(hΘk, sv.cv.fit$sieve.fit$beta_hat) %>% as.matrix
    }
    hΘ <- add_mats_pad(hΘ, hΘk[, -1] / 𝒟$K)
    estproxy  = lm(Yp ~ Xp + Zp - 1)
    hβp <- hβp + coef(estproxy) / 𝒟$K
  }
  M <- nrow(hΘ)
  Yt = 𝒟$targetdata$Y
  Xt = 𝒟$targetdata$X
  Ytest = 𝒟$testdata$Y
  Xtest = 𝒟$testdata$D$X
  
  #### Stage 1: Imputation - using the information of X
  t.sv.prp <- sieve_preprocess(Xt, M, type = 'cosine')
  tZt = t.sv.prp$Phi %*% hΘ
  
  #### Stage 2: Local debiasing
  tYt = c(Yt - Xt %*% hβp[1:p1] - tZt %*% hβp[(p1 + 1):p])
  combinedlassonew = cv.glmnet(cbind(Xt, tZt), tYt, alpha = 1, intercept = F)
  
  # Test errors
  test.sv.prp <- sieve_preprocess(Xtest, M, type = 'cosine')
  tZtest = test.sv.prp$Phi %*% hΘ
  Ypred = predict(combinedlassonew, newx = cbind(Xtest, tZtest)) +
    c(Xtest %*% hβp[1:p1] + tZtest %*% hβp[(p1 + 1):p])
  
  hdelta = coef(combinedlassonew)
  ee = norm(abs(𝒟$βt[1:p] - hdelta[1:p] - hβp[1:p]), '2')
  pe = norm(abs(Ytest - Ypred), '2')
  return(c(pe = pe, ee = ee))
}


# oracle ------------------------------------------------------------------


oracle = function(𝒟, Trgt) {
  p = 𝒟$p1 + 𝒟$p2
  Ytest = 𝒟$testdata$Y
  Xtest = 𝒟$testdata$D
  Ypred = as.matrix(Xtest) %*% Trgt$βt[1:p]
  pe = norm(Ytest - Ypred, '2')
  return(c(pe = pe, ee = 0))
}