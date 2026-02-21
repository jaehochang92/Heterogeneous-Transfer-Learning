# Codes by Max Russo (russo.325@osu.edu)
# Updated by Jae Chang (chang.2090@osu.edu)

library(curatedOvarianData)
library(glmnet)
library(dplyr)
library(limma)
library(ggplot2)
source('sv-cv.R')


# custom functions --------------------------------------------------------


calc_rmse = function(y_trth, y_prd) sqrt(mean((y_trth - y_prd)^2))
safe_ols_predict = function(x_trn, y_trn, x_tst) {
  x_trn_i = cbind(1, x_trn)
  x_tst_i = cbind(1, x_tst)
  qrx = qr(x_trn_i)
  keep = qrx$pivot[seq_len(qrx$rank)]
  coef_keep = coef(lm.fit(x = x_trn_i[, keep, drop = FALSE], y = y_trn))
  as.numeric(x_tst_i[, keep, drop = FALSE] %*% coef_keep)
}

# data prep ---------------------------------------------------------------


data(GSE26712_eset)
data(GSE9891_eset)
prxy_eset <- GSE26712_eset
trgt_eset <- GSE9891_eset

prxy_m = t(limma::as.matrix.ExpressionSet(prxy_eset))
trgt_m = t(limma::as.matrix.ExpressionSet(trgt_eset))

## select subjects id
id_sel_prxy = which(prxy_eset@phenoData$vital_status  =='deceased')
id_sel_trgt = which(trgt_eset@phenoData$vital_status =='deceased')

prxy_variables = read.csv("CuratedOvarianCancer/proxy_var.csv")[, 2] %>% sort
trgt_variables = read.csv("CuratedOvarianCancer/target_var.csv")[, 2] %>% sort


x_trgt = trgt_m[id_sel_trgt, trgt_variables]
y_trgt = log(trgt_eset@phenoData$days_to_death[id_sel_trgt])

x_prxy = prxy_m[id_sel_prxy, prxy_variables]
y_prxy = log(prxy_eset@phenoData$days_to_death[id_sel_prxy])

dim(x_prxy)
dim(x_trgt)

common_vars <- which(prxy_variables %in% trgt_variables)
extra_prxy_vars <- which(!prxy_variables %in% trgt_variables)
common_vars <- which(colnames(x_prxy) %in% colnames(x_trgt))
extra_prxy_vars = setdiff(1:ncol(x_prxy), common_vars)
p1 <- length(common_vars)
p2 <- length(extra_prxy_vars)
x_prxy <- x_prxy[, c(common_vars, extra_prxy_vars)]
x_trgt <- x_trgt[, colnames(x_prxy)[1:p1]]

if (p1 == 0) stop("No common variables between target and proxy datasets.")
if (p2 == 0) stop("Proxy dataset has no proxy-only variables after alignment.")


# feature mapping ---------------------------------------------------------


fmap_linear = cv_fitH_from_prxy(x_prxy[, 1:p1], x_prxy[, 1:p2 + p1], 'linear', 5)
fmap_sv = cv_fitH_from_prxy(x_prxy[, 1:p1], x_prxy[, 1:p2 + p1], 'sieve', 5)
save.image('CuratedOvarianCancer/021726.rdata')
# load('CuratedOvarianCancer/021726.rdata')

hat_z_prxy_linear <- fmap_predict(fmap_linear, x_prxy[, 1:p1], 'linear')
hat_z_prxy_sieve <- fmap_predict(fmap_sv, x_prxy[, 1:p1])
calc_rmse(x_prxy[, 1:p2 + p1], hat_z_prxy_linear)
calc_rmse(x_prxy[, 1:p2 + p1], hat_z_prxy_sieve)
par(mfrow = c(1, 2))
{
  i = 1
  plot(hat_z_prxy_linear[, i], x_prxy[, p1 + i])
  abline(0, 1)
  plot(hat_z_prxy_sieve[, i], x_prxy[, p1 + i])
  abline(0, 1)
}

cat("Linear prediction variance:", var(hat_z_prxy_linear[, 1]), "\n")
cat("Sieve prediction variance:", var(hat_z_prxy_sieve[, 1]), "\n")

# proxy model -------------------------------------------------------------

prxy_fit_cv = cv.glmnet(x_prxy, y_prxy, alpha = 1)
prxy_fit_cv_hmtl = cv.glmnet(x_prxy[, 1:p1], y_prxy, alpha = 1)
par(mfrow = c(1, 2))
plot(prxy_fit_cv)
plot(prxy_fit_cv_hmtl)


# CV ----------------------------------------------------------------------


n_trgt = nrow(x_trgt)
n_folds = 5
no <- 1234
{set.seed(no)
fold_indices = index_spliter(1:n_trgt, n_folds)
rmse_ols = rmse_htl_sv = rmse_htl_lin = rmse_hmtl = rmse_lasso = numeric(n_folds)
for (fold_idx in seq_len(n_folds)) {
  tst_idx = fold_indices[[fold_idx]]
  x_trn = x_trgt[-tst_idx, ]
  y_trn = y_trgt[-tst_idx]
  x_tst = x_trgt[tst_idx, ]
  y_tst = y_trgt[tst_idx]

  ## Stage 1: Imputation using X1
  ## linear
  hat_z_trn_lin = fmap_predict(fmap_linear, x_trn, 'linear')
  hat_z_tst_lin = fmap_predict(fmap_linear, x_tst, 'linear')
  ## sieve
  hat_z_trn_sv <- fmap_predict(fmap_sv, x_trn)
  hat_z_tst_sv <- fmap_predict(fmap_sv, x_tst)
  
  ## Stage 2: Estimation
  ## HTL
  hat_d_trn <- cbind(x_trn, hat_z_trn_lin)
  y_rsd_trn_lin = y_trn - predict(prxy_fit_cv, hat_d_trn, s = 'lambda.1se')
  htl_cv_lin = cv.glmnet(hat_d_trn, y_rsd_trn_lin, alpha = 1)
  hat_d_tst = cbind(x_tst, hat_z_tst_lin)
  y_prd_htl_lin = predict(htl_cv_lin, newx = hat_d_tst, s = "lambda.1se") + 
    predict(prxy_fit_cv, hat_d_tst, s = 'lambda.1se')
  
  hat_d_trn_sv <- cbind(x_trn, hat_z_trn_sv)
  y_rsd_trn_sv = y_trn - predict(prxy_fit_cv, hat_d_trn_sv, s = 'lambda.1se')
  htl_cv_sv = cv.glmnet(hat_d_trn_sv, y_rsd_trn_sv, alpha = 1)
  hat_d_tst_sv = cbind(x_tst, hat_z_tst_sv)
  y_prd_htl_sv = predict(htl_cv_sv, newx = hat_d_tst_sv, s = "lambda.1se") + 
    predict(prxy_fit_cv, hat_d_tst_sv, s = 'lambda.1se')

  ## HmTL (hmtl)
  y_rsd_hmtl = y_trn - predict(prxy_fit_cv_hmtl, x_trn, s = 'lambda.1se')
  hmtl_cv = cv.glmnet(x_trn, y_rsd_hmtl, alpha = 1)
  y_prd_hmtl = predict(hmtl_cv, x_tst, s = "lambda.1se") + 
    predict(prxy_fit_cv_hmtl, x_tst, s = "lambda.1se")
  
  ## lasso prediction
  lasso_cv = cv.glmnet(x_trn, y_trn, alpha = 1)
  y_prd_lasso = predict(lasso_cv, newx = x_tst, s = "lambda.1se")

  ## OLS with rank-deficiency handling
  y_prd_ols = safe_ols_predict(x_trn, y_trn, x_tst)

  ## errors
  rmse_ols[fold_idx] = calc_rmse(y_tst, y_prd_ols)
  rmse_htl_sv[fold_idx] = calc_rmse(y_tst, y_prd_htl_sv)
  rmse_htl_lin[fold_idx] = calc_rmse(y_tst, y_prd_htl_lin)
  rmse_hmtl[fold_idx] = calc_rmse(y_tst, y_prd_hmtl)
  rmse_lasso[fold_idx] = calc_rmse(y_tst, y_prd_lasso)
}

rmse_df = rbind(
  data.frame(fold = seq_len(n_folds), method = "HTL-sieve", rmse = rmse_htl_sv),
  data.frame(fold = seq_len(n_folds), method = "HTL-linear", rmse = rmse_htl_lin),
  data.frame(fold = seq_len(n_folds), method = "HmTL", rmse = rmse_hmtl),
  data.frame(fold = seq_len(n_folds), method = "Lasso", rmse = rmse_lasso),
  data.frame(fold = seq_len(n_folds), method = "OLS", rmse = rmse_ols)
) %>% filter(method != 'OLS')
rmse_df$method = factor(rmse_df$method, levels = c('HTL-sieve', "HTL-linear", "HmTL", "Lasso", "OLS"))

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
    inherit.aes = FALSE, shape = 95, size = 10, color = "black"
  ) +
  ggplot2::labs(x = NULL, y = "Root Mean Squared Error") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9)) +
  theme(axis.text.x = element_text(angle = 10, hjust = 1))

print(p)
ggplot2::ggsave("CuratedOvarianCancer/cv_results.pdf", plot = p, width = 4, height = 3.5)}