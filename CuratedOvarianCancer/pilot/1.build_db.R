library(dplyr)
library(tidyr)
library(SIS)

setwd('~/git/HTL_JASA_RnR/CuratedOvarianCancer/')
rm(list = ls())

ovarian_datasets <- readRDS('pilot/ovarian_datasets.rds')
data_pool <- lapply(ovarian_datasets, \(x) dim(x$X))
data_pool <- data_pool %>% do.call(rbind, .)
data_pool <- tibble(eset = rownames(data_pool), n = data_pool[, 1], d = data_pool[, 2])

# filtering
# data_pool <- data_pool |> filter(n > 40, d > 5e3)
# target.eset <- data_pool |> arrange(n, d) %>% .[1, 'eset', drop = T]
target.eset <- 'GSE9891'
proxy.eset <- data_pool$eset[data_pool$eset != target.eset]

target.eset <- ovarian_datasets[[target.eset]]
proxy.eset <- ovarian_datasets[proxy.eset]

# drop cols with NAs
lapply(proxy.eset, \(x) c(sum(is.na(x$X)), sum(is.na(x$y))))
drop_cols <- proxy.eset$GSE51088$X |> as_tibble() |>
  select(where(~any(is.na(.)))) |> names()
target.eset$X <- target.eset$X[, !colnames(target.eset$X) %in% drop_cols]

# target feature screening via SIS
sis <- SIS(target.eset$X, target.eset$y, 'gaussian', nsis = 3e3)
target.eset$X <- target.eset$X[, colnames(target.eset$X)[sis$sis.ix0]]

# search common genes
commons <- lapply(proxy.eset, \(x) intersect(colnames(target.eset$X), colnames(x$X)))
sapply(commons, length)
p1_genes <- colnames(target.eset$X)
for (d in proxy.eset) {
  p1_genes <- intersect(p1_genes, colnames(d$X))
}

# identify Z
unique_genes <- lapply(proxy.eset, \(x) setdiff(colnames(x$X), p1_genes))
unique_genes <- lapply(commons, \(x) setdiff(x, p1_genes))
set.seed(3)
unique_genes <- sample(unique_genes)
p2_genes <- colnames(proxy.eset$GSE13876$X)
for (i in seq_along(unique_genes)) {
  eset <- names(unique_genes)[i]
  new <- intersect(p2_genes, unique_genes[[i]])
  message(length(new))
  if (length(new) == 0) {
    proxy.eset[eset] <- NULL
    next
  }
  p2_genes <- new
}
p2_genes <- p2_genes[!p2_genes %in% drop_cols]
target.eset$X <- target.eset$X[, p1_genes]
p_genes <- c(p1_genes, p2_genes)
proxy.eset <- lapply(proxy.eset, \(x) {
  x$D$Z <- x$X[, p2_genes]
  x$D$X <- x$X[, p1_genes]
  x[c('D','y')]
  })
saveRDS(proxy.eset, 'pilot/proxy.rds')
saveRDS(target.eset, 'pilot/target.rds')