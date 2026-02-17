library(curatedOvarianData)
library(limma) ## needed for matrix manipulation
library(glmnet)


##  load targer and proxy data (98 target 26 proxy)
data(GSE9891_eset)
data(GSE26712_eset)

target_eset = t(as.matrix(GSE9891_eset))
proxy_eset = t(as.matrix(GSE26712_eset))

## select subjects id
id_sel_target = which(GSE9891_eset@phenoData$vital_status =='deceased')
id_sel_proxy = which(GSE26712_eset@phenoData$vital_status  =='deceased')

target_variables = read.csv("results/target_var.csv")[,2]
proxy_variables  = read.csv("results/proxy_var.csv")[,2]

X_target = target_eset[id_sel_target, which(colnames(target_eset) %in% target_variables)]
y_target = log(GSE9891_eset@phenoData$days_to_death[id_sel_target])

X_proxy = proxy_eset[id_sel_proxy, which(colnames(proxy_eset) %in% proxy_variables)]
y_proxy = log(GSE26712_eset@phenoData$days_to_death[id_sel_proxy])

## order proxy varible so that the first 69 variables are the commons one
ord_var = which(colnames(X_proxy) %in% colnames(X_target))
X_proxy = X_proxy[,c(ord_var,setdiff(1:NCOL(X_proxy),ord_var))]

## rergression on proxy data (one time)
## we use ridge regression for this step since p>n
proxy_ridge_cv = cv.glmnet(X_proxy, y_proxy, alpha = 0, intercept = FALSE)
proxy_beta = coef(proxy_ridge_cv, s = proxy_ridge_cv$lambda.min)[-1, 1]

proxyX2X1 = coef(lm(X_proxy[,-c(1:69)]~X_proxy[,1:69] -1))

## BASTANI ESTIMATOR
proxy_ridge_cv_bastani = cv.glmnet(X_proxy[,1:69],y_proxy, alpha = 0,intercept=FALSE)
proxy_beta_bastani = coef(proxy_ridge_cv_bastani, s = proxy_ridge_cv_bastani$lambda.min)[-1,1]


## create 200 bootstrap datasets and estimate coefficients 
set.seed(2175)
OLS = LASSO = BASTANI = TL = TL.sv = matrix(NA, 200, 69)
emap <- readRDS('results/emap.rds')
M <- nrow(emap)
for(b in 1:200) {
  b_id = sample(1:NROW(X_target), replace = TRUE)
  
  ## LASSO
  lasso_cv = cv.glmnet(X_target[b_id, ], y_target[b_id])
  lasso = coef(lasso_cv, s = 'lambda.min')[-1, 1]
  
  ## Bastani
  b_newYt = y_target[b_id] -  X_target[b_id, ] %*% proxy_beta_bastani
  blasso = cv.glmnet(X_target[b_id, ], b_newYt, alpha = 1)
  hdelta = coef(blasso, s = 'lambda.min')[, 1]
  bastani = hdelta[-1] + proxy_beta_bastani
  
  ## HTL w/ linear map
  b_X_imputed = X_target[b_id, ]  %*% proxyX2X1
  b_newYt = y_target[b_id] -  cbind(X_target[b_id, ], b_X_imputed) %*%  proxy_beta
  
  tl_lasso = cv.glmnet(cbind(X_target[b_id, ], b_X_imputed), y_target[b_id], alpha = 1)
  tl_hdelta = coef(tl_lasso, s = 'lambda.min')[, 1]
  
  TL_coef = tl_hdelta[-1] + proxy_beta
  TL_coef = TL_coef[1:69]
  
  ## HTL w/ sieve
  sv_target <- sieve_preprocess(X_target[b_id, ], M)
  b_X_imputed = sv_target$Phi %*% emap
  b_newYt = y_target[b_id] - cbind(X_target[b_id, ], b_X_imputed) %*% proxy_beta
  
  tl_lasso = cv.glmnet(cbind(X_target[b_id, ], b_X_imputed), y_target[b_id], alpha = 1)
  tl_hdelta = coef(tl_lasso, s = 'lambda.min')[, 1]
  
  TL.sv_coef = tl_hdelta[-1] + proxy_beta
  TL.sv_coef = TL.sv_coef[1:69]
  
  ## save results
  OLS[b, ]  = ols
  LASSO[b, ] = lasso
  BASTANI[b, ] = bastani
  TL[b, ] = TL_coef
  TL.sv[b, ] = TL.sv_coef
}

pdf("bootstrap.pdf")
mar.default <- c(5, 4, 4, 2)
par(mfrow = c(2, 3), mar =  mar.default + c(3, 0, 0, 0))
for (d in 1:6) {
  boxplot(
    TL.sv[, d],
    TL[, d],
    BASTANI[, d],
    LASSO[, d],
    names = c('Htg.TL.sieve', 'Htg.TL.linear', 'Hmg.TL', 'Lasso'),
    main = colnames(X_target)[d],
    cex.lab = 1.2,
    boxwex = .3,
    xaxt = 'n'
  )
  axis(1,
       1:4,
       c('Htg.TL.sieve', 'Htg.TL.linear', 'Hmg.TL', 'Lasso'),
       las = 2)
}
dev.off()