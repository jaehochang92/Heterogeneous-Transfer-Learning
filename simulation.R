## --- Single source of truth: vectors + derived constants ---


Kv  <- c(2, 4, 8)[]
ntv <- c(30, 1e2, 2e2)[]
p1v <- c(12, 25, 50)[]
p2v <- c(12, 25, 50)[]
pv  <- p1v + p2v
npv <- 2 * pv
np0 <- max(npv); nt0 <- min(ntv); p10 <- min(p1v); p20 <- min(p2v)

Nrep  <- 2e2
ntest <- 2e2

# regime ------------------------------------------------------------------


regime <- function(seedno, K, p1, p2) {
  p <- p1 + p2
  S <- sqrt(p / 2) %>% round
  Prxy <- list(
    Σ1p = list(),
    Σξp = list(),
    δst_sprs = list(),
    δst_nsprs = list()
  )
  Trgt <- list()
  set.seed(seedno)
  for (k in 1:K) {
    Prxy$Σ1p[[k]] <- diag(p1)
    Prxy$Σξp[[k]] <- diag(p2)
    B1 <- c(rep(1, S), rep(0, p1 - S)) %>% sample
    B2 <- c(rep(1, S), rep(0, p2 - S)) %>% sample
    contrast <- rnorm(p, 0, 1 / 4)
    Prxy$δst_sprs[[k]] <- contrast * c(B1, B2) # sparse case
    Prxy$δst_nsprs[[k]] <- contrast # non-sparse case
  }
  b <- 15
  Trgt$βt <- rep(c(1, rep(0, p / b)), b)[1:p]
  Trgt$Σξt <- diag(p2)
  return(list(Prxy = Prxy, Trgt = Trgt))
}

Kmax = max(Kv) 
p1max = max(p1v) 
p2max = max(p2v)
regime0 <- regime(1595, Kmax, p1max, p2max) # initialize regime

# single run --------------------------------------------------------------


single_sim <- function(K, np, nt, p1, p2, ntest, Prxy, Trgt) {
  𝒟 <- DGP_transfer(K, np, nt, p1, p2, ntest, Prxy, Trgt)
  output <- tibble()
  for (δsti in c('sp', 'nsp')) {
    output <- bind_rows(
      output,
      tibble(
        lasso_trgt = lasso_trgt(𝒟),
        TL_matched = TL_matched(𝒟, δsti),
        TL_impute = TL_impute(𝒟, δsti),
        oracle = oracle(𝒟, Trgt),
        loss = c('PrdErr', 'EstErr'),
        sparse = δsti
      )
    )
  }
  return(output)
}

# Main simulation ---------------------------------------------------------


main <- function(rep, K, np, nt, p1, p2, ntest, regime0) {
  # Monte Carlo
  mcrep <- pbreplicate(rep,
                       single_sim(K, np, nt, p1, p2, ntest,
                                  regime0$Prxy, regime0$Trgt),
                       F)
  mcrep <- do.call("bind_rows", mcrep)
  # aggregation and ggplot
  gg_arr <- mcrep %>% pivot_longer(-c(sparse, loss),
                                   names_to = 'Methodology',
                                   values_to = 'value') %>%
    mutate(
      dlt_sprs = factor(
        sparse,
        levels = c('sp', 'nsp'),
        labels = c('Y', 'N')
      ),
      Methodology = factor(
        Methodology,
        levels = c('oracle', 'TL_impute', 'TL_matched', 'lasso_trgt')
      )
    )
  return(gg_arr)
}