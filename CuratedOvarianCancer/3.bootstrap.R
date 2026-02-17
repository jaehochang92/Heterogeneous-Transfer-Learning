library(curatedOvarianData)
library(glmnet)

## load target and proxy data
data(GSE9891_eset)
data(GSE26712_eset)
target_eset = t(as.matrix(GSE9891_eset))
proxy_eset = t(as.matrix(GSE26712_eset))

## select subject ids
target_deceased_ids = which(GSE9891_eset@phenoData$vital_status == "deceased")
proxy_deceased_ids = which(GSE26712_eset@phenoData$vital_status == "deceased")

target_variables = read.csv("results/target_var.csv")[, 2]
proxy_variables = read.csv("results/proxy_var.csv")[, 2]

x_target = target_eset[target_deceased_ids, which(colnames(target_eset) %in% target_variables)]
y_target = log(GSE9891_eset@phenoData$days_to_death[target_deceased_ids])
x_proxy = proxy_eset[proxy_deceased_ids, which(colnames(proxy_eset) %in% proxy_variables)]
y_proxy = log(GSE26712_eset@phenoData$days_to_death[proxy_deceased_ids])

## align shared and proxy-only variables
common_proxy_idx = which(colnames(x_proxy) %in% colnames(x_target))
x_proxy = x_proxy[, c(common_proxy_idx, setdiff(seq_len(ncol(x_proxy)), common_proxy_idx))]
n_common = length(common_proxy_idx)

## regression on proxy data (one time)
proxy_ridge_cv = cv.glmnet(x_proxy, y_proxy, alpha = 0, intercept = FALSE)
proxy_beta = coef(proxy_ridge_cv, s = proxy_ridge_cv$lambda.min)[-1, 1]
proxy_x2_on_x1_coef = coef(lm(x_proxy[, -(1:n_common), drop = FALSE] ~ x_proxy[, 1:n_common, drop = FALSE] - 1))

## BASTANI estimator
proxy_ridge_cv_bastani = cv.glmnet(x_proxy[, 1:n_common, drop = FALSE], y_proxy, alpha = 0, intercept = FALSE)
proxy_beta_bastani = coef(proxy_ridge_cv_bastani, s = proxy_ridge_cv_bastani$lambda.min)[-1, 1]

## create bootstrap datasets and estimate coefficients
set.seed(2175)
n_boot = 200
lasso_coef_mat = bastani_coef_mat = htl_linear_coef_mat = htl_sieve_coef_mat = matrix(NA_real_, n_boot, n_common)
emap = readRDS("results/emap.rds")
sieve_basis_n = nrow(emap)

for (boot_iter in seq_len(n_boot)) {
  boot_ids = sample(seq_len(nrow(x_target)), replace = TRUE)
  x_boot = x_target[boot_ids, , drop = FALSE]
  y_boot = y_target[boot_ids]

  ## LASSO
  lasso_cv = cv.glmnet(x_boot, y_boot, alpha = 1)
  lasso_coef = coef(lasso_cv, s = "lambda.min")[-1, 1]

  ## Bastani
  new_y_bastani = y_boot - x_boot %*% proxy_beta_bastani
  bastani_lasso = cv.glmnet(x_boot, new_y_bastani, alpha = 1)
  delta_bastani = coef(bastani_lasso, s = "lambda.min")[-1, 1]
  bastani_coef = delta_bastani + proxy_beta_bastani

  ## HTL with linear map
  x_imputed_linear = x_boot %*% proxy_x2_on_x1_coef
  new_y_linear = y_boot - cbind(x_boot, x_imputed_linear) %*% proxy_beta
  htl_linear_lasso = cv.glmnet(cbind(x_boot, x_imputed_linear), new_y_linear, alpha = 1)
  htl_linear_coef = (coef(htl_linear_lasso, s = "lambda.min")[-1, 1] + proxy_beta)[1:n_common]

  ## HTL with sieve
  sieve_target = sieve_preprocess(x_boot, sieve_basis_n)
  x_imputed_sieve = sieve_target$Phi %*% emap
  new_y_sieve = y_boot - cbind(x_boot, x_imputed_sieve) %*% proxy_beta
  htl_sieve_lasso = cv.glmnet(cbind(x_boot, x_imputed_sieve), new_y_sieve, alpha = 1)
  htl_sieve_coef = (coef(htl_sieve_lasso, s = "lambda.min")[-1, 1] + proxy_beta)[1:n_common]

  ## save results
  lasso_coef_mat[boot_iter, ] = lasso_coef
  bastani_coef_mat[boot_iter, ] = bastani_coef
  htl_linear_coef_mat[boot_iter, ] = htl_linear_coef
  htl_sieve_coef_mat[boot_iter, ] = htl_sieve_coef
}

pdf("bootstrap.pdf")
mar_default = c(5, 4, 4, 2)
par(mfrow = c(2, 3), mar = mar_default + c(3, 0, 0, 0))
for (plot_idx in 1:6) {
  boxplot(
    htl_sieve_coef_mat[, plot_idx],
    htl_linear_coef_mat[, plot_idx],
    bastani_coef_mat[, plot_idx],
    lasso_coef_mat[, plot_idx],
    names = c("Htg.TL.sieve", "Htg.TL.linear", "Hmg.TL", "Lasso"),
    main = colnames(x_target)[plot_idx],
    cex.lab = 1.2,
    boxwex = 0.3,
    xaxt = "n"
  )
  axis(1, 1:4, c("Htg.TL.sieve", "Htg.TL.linear", "Hmg.TL", "Lasso"), las = 2)
}
dev.off()
