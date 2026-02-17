library(curatedOvarianData)
library(limma) ## needed for matrix manipulation
library(glmnet)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Please install it with install.packages('ggplot2').")
}

## load target and proxy data
data(GSE9891_eset)
data(GSE26712_eset)
target_eset = t(as.matrix(GSE9891_eset))
proxy_eset = t(as.matrix(GSE26712_eset))

## select subjects id
id_sel_target = which(GSE9891_eset@phenoData$vital_status == "deceased")
id_sel_proxy = which(GSE26712_eset@phenoData$vital_status == "deceased")

read_var_file = function(filename) {
  if (file.exists(filename)) return(read.csv(filename)[, 2])
  read.csv(file.path("CuratedOvarianCancer", filename))[, 2]
}
target_variables = read_var_file("target_26_var.csv")
proxy_variables = read_var_file("proxy_98_var.csv")

X_target = target_eset[id_sel_target, which(colnames(target_eset) %in% target_variables)]
y_target = log(GSE9891_eset@phenoData$days_to_death[id_sel_target])
X_proxy = proxy_eset[id_sel_proxy, which(colnames(proxy_eset) %in% proxy_variables)]
y_proxy = log(GSE26712_eset@phenoData$days_to_death[id_sel_proxy])

## align common variables by name for both datasets
common_vars = intersect(colnames(X_target), colnames(X_proxy))
extra_proxy_vars = setdiff(colnames(X_proxy), common_vars)
X_target = scale(X_target[, common_vars, drop = FALSE], scale = F)
X_proxy = scale(X_proxy[, c(common_vars, extra_proxy_vars), drop = FALSE], scale = F)

n_common = length(common_vars)
if (n_common == 0) stop("No common variables between target and proxy datasets.")
if (ncol(X_proxy) <= n_common) stop("Proxy dataset has no proxy-only variables after alignment.")

## cross validation prediction error
set.seed(13)
n_folds = 5
n_target = nrow(X_target)
fold_id = sample(rep(1:n_folds, length.out = n_target))
fold_indices = split(seq_len(n_target), fold_id)
rmse = function(y_true, y_pred) sqrt(mean((y_true - y_pred)^2))

safe_ols_predict = function(X_train, y_train, X_test) {
  X_train_i = cbind(1, X_train)
  X_test_i = cbind(1, X_test)
  qrx = qr(X_train_i)
  keep = qrx$pivot[seq_len(qrx$rank)]
  coef_keep = coef(lm.fit(x = X_train_i[, keep, drop = FALSE], y = y_train))
  as.numeric(X_test_i[, keep, drop = FALSE] %*% coef_keep)
}

## results place holders
RMSE_t = RMSE_im = RMSE_lasso = RMSE_b = numeric(n_folds)

## regression on proxy data (one time)
proxy_ridge_cv = cv.glmnet(X_proxy, y_proxy, alpha = 0, intercept = FALSE)
proxy_beta = coef(proxy_ridge_cv, s = proxy_ridge_cv$lambda.min)[-1, 1]
proxyX2X1 = coef(lm(X_proxy[, -(1:n_common), drop = FALSE] ~ X_proxy[, 1:n_common, drop = FALSE] - 1))
proxy_ridge_cv_bastani = cv.glmnet(X_proxy[, 1:n_common, drop = FALSE], y_proxy, alpha = 0, intercept = FALSE)
proxy_beta_bastani = coef(proxy_ridge_cv_bastani, s = proxy_ridge_cv_bastani$lambda.min)[-1, 1]

## CV
for (k in seq_len(n_folds)) {
  test = fold_indices[[k]]
  train = setdiff(1:n_target, test)
  cv_Xtarget = X_target[train, ]
  cv_ytarget = y_target[train]
  X_test = X_target[test, , drop = FALSE]
  y_test = y_target[test]

  ## Stage 1: Imputation using X1
  cv_X_imputed = cv_Xtarget %*% proxyX2X1
  cv_X_imputed_test = X_test %*% proxyX2X1

  ## Stage 2: Estimation
  cv_newY_target = cv_ytarget - cbind(cv_Xtarget, cv_X_imputed) %*% proxy_beta
  cv_combinedlassonew = cv.glmnet(cbind(cv_Xtarget, cv_X_imputed), cv_newY_target, alpha = 1)
  cv_newX = cbind(X_test, cv_X_imputed_test)
  TL_Ypred = predict(cv_combinedlassonew, newx = cv_newX, s = "lambda.min") + cv_newX %*% proxy_beta

  ## lasso prediction
  lcv = cv.glmnet(cv_Xtarget, cv_ytarget, alpha = 1)
  lasso_Ypred = predict(lcv, newx = X_test, s = "lambda.min")

  ## OLS with rank-deficiency handling
  ols_Ypred = safe_ols_predict(cv_Xtarget, cv_newY_target, X_test)

  ## BASTANI
  newYt = cv_ytarget - cv_Xtarget %*% proxy_beta_bastani
  blasso = cv.glmnet(cv_Xtarget, newYt, alpha = 1)
  B_Ypred = predict(blasso, newx = X_test, s = "lambda.min") + X_test %*% proxy_beta_bastani

  ## errors
  RMSE_t[k] = rmse(y_test, ols_Ypred)
  RMSE_im[k] = rmse(y_test, TL_Ypred)
  RMSE_lasso[k] = rmse(y_test, lasso_Ypred)
  RMSE_b[k] = rmse(y_test, B_Ypred)
}

rmse_df = rbind(
  data.frame(fold = seq_len(n_folds), method = "HTL", rmse = RMSE_im),
  data.frame(fold = seq_len(n_folds), method = "HmTL", rmse = RMSE_b),
  data.frame(fold = seq_len(n_folds), method = "Lasso", rmse = RMSE_lasso),
  data.frame(fold = seq_len(n_folds), method = "OLS", rmse = RMSE_t)
)
rmse_df$method = factor(rmse_df$method, levels = c("HTL", "HmTL", "Lasso", "OLS"))

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
  ggplot2::geom_point(
    shape = 4, size = 1.2, stroke = 1, color = "red", alpha = 0.9
    # , position = ggplot2::position_jitter(width = 0.04, height = 0)
  ) +
  ggplot2::geom_errorbar(
    data = rmse_summary,
    ggplot2::aes(x = method, ymin = lower, ymax = upper),
    inherit.aes = FALSE, width = 0, linewidth = 0.5, color = "black"
  ) +
  ggplot2::geom_point(
    data = rmse_summary,
    ggplot2::aes(x = method, y = mean),
    inherit.aes = FALSE, shape = 95, size = 1, color = "black"
  ) +
  ggplot2::labs(x = NULL, y = "Root Mean Squared Error") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9))

print(p)
ggplot2::ggsave("cv_results.pdf", plot = p, width = 3, height = 2.5)
