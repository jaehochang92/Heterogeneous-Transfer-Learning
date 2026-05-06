# Codes by Max Russo (russo.325@osu.edu)
# Updated by Jae Chang (chang.2090@osu.edu)

library(curatedOvarianData)
library(glmnet)
library(dplyr)
library(limma)
library(ggplot2)


# helper functions -------------------------------------------------------

calc_rmse = function(y_trth, y_prd) sqrt(mean((y_trth - y_prd)^2))

safe_ols_predict = function(x_trn, y_trn, x_tst) {
  x_trn_i = cbind(1, x_trn)
  x_tst_i = cbind(1, x_tst)
  qrx = qr(x_trn_i)
  keep = qrx$pivot[seq_len(qrx$rank)]
  coef_keep = coef(lm.fit(x = x_trn_i[, keep, drop = FALSE], y = y_trn))
  as.numeric(x_tst_i[, keep, drop = FALSE] %*% coef_keep)
}

fit_htl <- function(x_trn, y_trn, x_tst, hat_z_trn = NULL, hat_z_tst = NULL, baseline_model) {
  hat_d_trn <- if (is.null(hat_z_trn)) x_trn else cbind(x_trn, hat_z_trn)
  hat_d_tst <- if (is.null(hat_z_tst)) x_tst else cbind(x_tst, hat_z_tst)
  baseline_pred_trn <- predict(baseline_model, newx = hat_d_trn, s = 'lambda.min')
  baseline_pred_tst <- predict(baseline_model, newx = hat_d_tst, s = 'lambda.min')
  htl_cv <- cv.glmnet(x = hat_d_trn, y = y_trn, offset = baseline_pred_trn, alpha = 1)
  predict(htl_cv, newx = hat_d_tst, newoffset = baseline_pred_tst, s = 'lambda.min')
}

# data prep ---------------------------------------------------------------


data(GSE26712_eset)
data(GSE9891_eset)
prxy_eset <- GSE26712_eset
trgt_eset <- GSE9891_eset
# prxy_eset <- GSE9891_eset
# trgt_eset <- GSE26712_eset

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
np <- nrow(x_prxy)
p1 <- length(common_vars)
p2 <- length(extra_prxy_vars)
x_prxy <- x_prxy[, c(common_vars, extra_prxy_vars)]
x_trgt <- x_trgt[, colnames(x_prxy)[1:p1]]

if (p1 == 0) stop("No common variables between target and proxy datasets.")
if (p2 == 0) stop("Proxy dataset has no proxy-only variables after alignment.")

prxy_min <- apply(x_prxy[, 1:p1], 2, min)
prxy_max <- apply(x_prxy[, 1:p1], 2, max)
fmap_predict <- function(fmap, X, model = 'sieve', base_min = prxy_min, base_max = prxy_max) {
  if (model == 'linear') return(cbind(1, X) %*% fmap)
  X_decoy <- rbind(base_min, base_max, X)
  X_sv_prp <- sieve_preprocess(X_decoy, nrow(fmap), type = 'cosine')
  Phi_decoy <- X_sv_prp$Phi
  X_sv_prp <- sieve_preprocess(X, nrow(fmap), type = 'cosine')
  Phi <- X_sv_prp$Phi
  Phi %*% fmap
}

# feature mapping ---------------------------------------------------------


source('sv-cv.R')
fmap_linear = cv_fitH_from_prxy(x_prxy[, 1:p1], x_prxy[, 1:p2 + p1], 'linear')
fmap_sv = cv_fitH_from_prxy(x_prxy[, 1:p1], x_prxy[, 1:p2 + p1], 'sieve')
# save.image('CuratedOvarianCancer/fmapped.rdata')
load('CuratedOvarianCancer/fmapped.rdata')

hat_z_prxy_linear <- fmap_predict(fmap_linear, x_prxy[, 1:p1], 'linear')
hat_z_prxy_sieve <- fmap_predict(fmap_sv, x_prxy[, 1:p1])
calc_rmse(x_prxy[, 1:p2 + p1], hat_z_prxy_linear)
calc_rmse(x_prxy[, 1:p2 + p1], hat_z_prxy_sieve)

cat("Linear prediction variance:", var(hat_z_prxy_linear[, 1]), "\n")
cat("Sieve prediction variance:", var(hat_z_prxy_sieve[, 1]), "\n")

# proxy model -------------------------------------------------------------

elnet_alp <- 0
prxy_fit_cv <- cv.glmnet(x_prxy, y_prxy, alpha = elnet_alp)
prxy_fit_cv_hmtl <- cv.glmnet(x_prxy[, 1:p1], y_prxy, alpha = elnet_alp)
par(mfrow = c(1, 2))
plot(prxy_fit_cv, main = 'Full Model')
plot(prxy_fit_cv_hmtl)
htl_cf <- coef(prxy_fit_cv, 'lambda.min')
hmtl_cf <- coef(prxy_fit_cv_hmtl, 'lambda.min')
padder0(HTL = htl_cf, HmTL = hmtl_cf)[1:p1 + 1, ]

# CV ----------------------------------------------------------------------


n_trgt = nrow(x_trgt)
n_folds = 5
run.once(2175)
run.once <- function(no) {
  set.seed(no)
  fold_indices = index_spliter(1:n_trgt, n_folds)
  
  results <- lapply(seq_len(n_folds), function(fold_idx) {
    tst_idx = fold_indices[[fold_idx]]
    x_trn = x_trgt[-tst_idx, ]; y_trn = y_trgt[-tst_idx]
    x_tst = x_trgt[tst_idx, ];  y_tst = y_trgt[tst_idx]
    
    ## Imputation
    hat_z_trn_lin = fmap_predict(fmap_linear, x_trn, 'linear')
    hat_z_tst_lin = fmap_predict(fmap_linear, x_tst, 'linear')
    hat_z_trn_sv = fmap_predict(fmap_sv, x_trn)
    hat_z_tst_sv = fmap_predict(fmap_sv, x_tst)
    
    ## Predictions
    predictions <- list(
      `HTL-sieve` = fit_htl(x_trn, y_trn, x_tst, hat_z_trn_sv, hat_z_tst_sv, prxy_fit_cv),
      `HTL-linear` = fit_htl(x_trn, y_trn, x_tst, hat_z_trn_lin, hat_z_tst_lin, prxy_fit_cv),
      HmTL = fit_htl(x_trn, y_trn, x_tst, baseline_model = prxy_fit_cv_hmtl),
      Lasso = predict(cv.glmnet(x_trn, y_trn, alpha = 1), newx = x_tst, s = "lambda.min")
    )
    
    data.frame(
      fold = fold_idx,
      method = names(predictions),
      rmse = sapply(predictions, calc_rmse, y_trth = y_tst),
      row.names = NULL
    )
  }) %>% bind_rows() %>%
    mutate(method = factor(method, levels = c('HTL-sieve', 'HTL-linear', 'HmTL', 'Lasso')))
  
  summary_stats <- results %>%
    group_by(method) %>%
    summarise(mean = mean(rmse), se = sd(rmse) / sqrt(n()), .groups = 'drop') %>%
    mutate(lower = mean - 2 * se, upper = mean + 2 * se)
  
  ggplot(results, aes(method, rmse)) +
    geom_hline(yintercept = pretty(range(results$rmse), 6), 
               color = "grey85", linewidth = 0.3) +
    geom_point(shape = 4, size = 1, stroke = 1, color = "red", alpha = 0.9) +
    geom_errorbar(aes(y = mean, ymin = lower, ymax = upper), 
                  summary_stats, width = 0, linewidth = 0.5) +
    geom_point(aes(y = mean), summary_stats, shape = 95, size = 10) +
    labs(x = NULL, y = "Root Mean Squared Error") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(size = 9, angle = 10, hjust = 1))
  
  ggsave("CuratedOvarianCancer/cv_results.pdf", width = 4, height = 3.5)
}