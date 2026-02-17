# Codes by Max Russo (russo.325@osu.edu)
# Updated by Jae Chang (chang.2090@osu.edu)

# preliminary functions
calc_rmse = function(y_trth, y_prd) sqrt(mean((y_trth - y_prd)^2))
safe_ols_predict = function(x_trn, y_trn, x_tst) {
  x_trn_i = cbind(1, x_trn)
  x_tst_i = cbind(1, x_tst)
  qrx = qr(x_trn_i)
  keep = qrx$pivot[seq_len(qrx$rank)]
  coef_keep = coef(lm.fit(x = x_trn_i[, keep, drop = FALSE], y = y_trn))
  as.numeric(x_tst_i[, keep, drop = FALSE] %*% coef_keep)
}

library(curatedOvarianData)
library(limma) ## needed for matrix manipulation
library(glmnet)
library(dplyr)
library(ggplot2)

source('sv-cv.R')
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Please install it with install.packages('ggplot2').")
}
## load target and proxy data
data(GSE9891_eset)
data(GSE26712_eset)
target_eset = t(as.matrix(GSE9891_eset))
prxyeset = t(as.matrix(GSE26712_eset))

## select subject ids
target_deceased_ids = which(GSE9891_eset@phenoData$vital_status == "deceased")
prxydeceased_ids = which(GSE26712_eset@phenoData$vital_status == "deceased")

read_var_file = function(filename) {
  if (file.exists(filename)) return(read.csv(filename)[, 2])
  read.csv(file.path("CuratedOvarianCancer", filename))[, 2]
}
target_variables = read_var_file("CuratedOvarianCancer/target_26_var.csv")
prxyvariables = read_var_file("CuratedOvarianCancer/proxy_98_var.csv")

x_trgt = target_eset[target_deceased_ids, which(colnames(target_eset) %in% target_variables)]
y_trgt = log(GSE9891_eset@phenoData$days_to_death[target_deceased_ids])
x_proxy = prxyeset[prxydeceased_ids, which(colnames(prxyeset) %in% prxyvariables)]
y_proxy = log(GSE26712_eset@phenoData$days_to_death[prxydeceased_ids])

## align common variables by name for both datasets
common_vars = intersect(colnames(x_trgt), colnames(x_proxy))
p1 <- length(common_vars)
extra_prxyvars = setdiff(colnames(x_proxy), common_vars)
p2 <- length(extra_prxyvars)
x_trgt = scale(x_trgt[, common_vars, drop = FALSE], scale = FALSE)
x_proxy = scale(x_proxy[, c(common_vars, extra_prxyvars), drop = FALSE], scale = FALSE)

if (p1 == 0) stop("No common variables between target and proxy datasets.")
if (p2 == 0) stop("Proxy dataset has no proxy-only variables after alignment.")


## regression on proxy data (one time)
prxyridge_cv = cv.glmnet(x_proxy, y_proxy, alpha = 0, intercept = T)
prxyridge_cv_hmtl = cv.glmnet(x_proxy[, 1:p1, drop = FALSE], y_proxy, alpha = 0, intercept = T)
prxybeta = coef(prxyridge_cv, s = prxyridge_cv$lambda.min)[-1, 1]
prxybeta_hmtl = coef(prxyridge_cv_hmtl, s = prxyridge_cv_hmtl$lambda.min)[-1, 1]

## feature map
prxyH_lin = cv_fitH_from_prxy(x_proxy, p1, 'linear')
prxyH_sv = cv_fitH_from_prxy(x_proxy, p1, 'sieve', 5, 8)

# save.image('CuratedOvarianCancer/021726.rdata')
load('CuratedOvarianCancer/021726.rdata')

seed.no <- sample(1:3e4, 1);{set.seed(seed.no)
n_folds = 5
n_trgt = nrow(x_trgt)
fold_id = sample(rep(1:n_folds, length.out = n_trgt))
fold_indices = split(seq_len(n_trgt), fold_id)

## CV
rmse_ols = rmse_htl_sv = rmse_htl_lin = rmse_hmtl = rmse_lasso = numeric(n_folds)
for (fold_idx in seq_len(n_folds)) {
  tst_idx = fold_indices[[fold_idx]]
  trn_idx = setdiff(1:n_trgt, tst_idx)
  x_trn = x_trgt[trn_idx, ]
  y_trn = y_trgt[trn_idx]
  x_tst = x_trgt[tst_idx, , drop = FALSE]
  y_tst = y_trgt[tst_idx]

  ## Stage 1: Imputation using X1
  hat_z_trn_lin = x_trn %*% prxyH_lin
  hat_z_tst_lin = x_tst %*% prxyH_lin
  M <- nrow(prxyH_sv)
  Phi <- sieve_preprocess(x_trn, basisN = M, type = 'cosine')$Phi
  hat_z_trn_sv = Phi %*% prxyH_sv
  Phi <- sieve_preprocess(x_tst, basisN = M, type = 'cosine')$Phi
  hat_z_tst_sv = Phi %*% prxyH_sv

  ## Stage 2: Estimation
  ## HTL
  hat_d_trn <- cbind(x_trn, hat_z_trn_lin)
  y_rsd_trn_lin = y_trn - hat_d_trn %*% prxybeta
  htl_cv_lin = cv.glmnet(hat_d_trn, y_rsd_trn_lin, alpha = 1, intercept = T)
  hat_d_test = cbind(x_tst, hat_z_tst_lin)
  y_prd_htl_lin = predict(htl_cv_lin, newx = hat_d_test, s = "lambda.min") + hat_d_test %*% prxybeta
  
  hat_d_tst <- cbind(x_tst, hat_z_tst_lin)
  y_rsd_trn_sv = y_trn - cbind(x_trn, hat_z_trn_sv) %*% prxybeta
  htl_cv_sv = cv.glmnet(cbind(x_trn, hat_z_trn_sv), y_rsd_trn_sv, alpha = 1)
  hat_d_test = cbind(x_tst, hat_z_tst_sv)
  y_prd_htl_sv = predict(htl_cv_sv, newx = hat_d_test, s = "lambda.min") + hat_d_test %*% prxybeta

  ## lasso prediction
  lasso_cv = cv.glmnet(x_trn, y_trn, alpha = 1, intercept = T)
  y_prd_lasso = predict(lasso_cv, newx = x_tst, s = "lambda.min")

  ## OLS with rank-deficiency handling
  y_prd_ols = safe_ols_predict(x_trn, y_rsd_trn_lin, x_tst)

  ## HmTL (hmtl)
  y_residual_hmtl = y_trn - x_trn %*% prxybeta_hmtl
  hmtl_cv = cv.glmnet(x_trn, y_residual_hmtl, alpha = 1, intercept = T)
  y_prd_hmtl = predict(hmtl_cv, newx = x_tst, s = "lambda.min") + x_tst %*% prxybeta_hmtl

  ## errors
  rmse_ols[fold_idx] = calc_rmse(y_tst, y_prd_ols)
  rmse_htl_sv[fold_idx] = calc_rmse(y_tst, y_prd_htl_sv)
  rmse_htl_lin[fold_idx] = calc_rmse(y_tst, y_prd_htl_lin)
  rmse_lasso[fold_idx] = calc_rmse(y_tst, y_prd_lasso)
  rmse_hmtl[fold_idx] = calc_rmse(y_tst, y_prd_hmtl)
}

rmse_df = rbind(
  data.frame(fold = seq_len(n_folds), method = "HTL-sieve", rmse = rmse_htl_sv),
  data.frame(fold = seq_len(n_folds), method = "HTL-linear", rmse = rmse_htl_lin),
  data.frame(fold = seq_len(n_folds), method = "HmTL", rmse = rmse_hmtl),
  data.frame(fold = seq_len(n_folds), method = "Lasso", rmse = rmse_lasso),
  data.frame(fold = seq_len(n_folds), method = "OLS", rmse = rmse_ols)
)
rmse_df$method = factor(rmse_df$method, levels = c('HTL-sieve', "HTL-linear", "HmTL", "Lasso", "OLS"))
rmse_df <- rmse_df %>% filter(method != 'OLS')

rmse_summary = aggregate(rmse ~ method, rmse_df, function(x) {
  c(
    mean = mean(x),
    lower = mean(x) - 2 * sd(x) / sqrt(length(x)),
    upper = mean(x) + 2 * sd(x) / sqrt(length(x))
  )
})

rmse_summary = data.frame(
  method = rmse_summary$method,
  mean = rmse_summary$rmse[, "mean"],
  lower = rmse_summary$rmse[, "lower"],
  upper = rmse_summary$rmse[, "upper"]
)

p = ggplot2::ggplot(rmse_df, ggplot2::aes(x = method, y = rmse)) +
  ggplot2::geom_hline(yintercept = pretty(range(rmse_df$rmse), n = 6), color = "grey85", linewidth = 0.3) +
  ggplot2::geom_point(shape = 4, size = 1, stroke = 1, color = "red", alpha = 0.9) +
  ggplot2::geom_errorbar(
    data = rmse_summary,
    ggplot2::aes(x = method, ymin = lower, ymax = upper),
    inherit.aes = FALSE, width = 0, linewidth = 0.5, color = "black"
  ) +
  ggplot2::geom_point(
    data = rmse_summary,
    ggplot2::aes(x = method, y = mean),
    inherit.aes = FALSE, shape = 95, size = 4, color = "black"
  ) +
  ggplot2::labs(x = NULL, y = "Root Mean Squared Error") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9)) +
  theme(axis.text.x = element_text(angle = 10, hjust = 1))

print(p)
ggplot2::ggsave("CuratedOvarianCancer/cv_results.pdf", plot = p, width = 4, height = 3.5)}