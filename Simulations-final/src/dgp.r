# Configure data distributions --------------------------------------------

dgp.confg <- function(seedno, K, p1, p2, fmap) {
  set.seed(seedno)
  p <- p1 + p2
  s0 <- sqrt(p / 2) %>% round
  Prxy <- list(
    Σ1p = list(), Σξp = list(), δst_sprs = list(), δst_nsprs = list()
  )
  Trgt <- list()
  if (fmap == 'linear') {
    Trgt$Θt <- matrix(rbeta(p1 * p2, 10, 10), p1) * 10
  }
  for (k in 1:K) {
    Prxy$Σ1p[[k]] <- diag(p1)
    Prxy$Σξp[[k]] <- diag(p2)
    if (fmap == 'linear') {
      Prxy$Θp[[k]] <- Trgt$Θt + 1 / 3 * (matrix(rbeta(p1 * p2, 4, 4), p1) - 1 / 2)
    }
    β1 <- c(rep(1, s0), rep(0, p1 - s0)) %>% sample
    β2 <- c(rep(1, s0), rep(0, p2 - s0)) %>% sample
    contrast <- rnorm(p, 0, 1 / 4)
    Prxy$δst_sprs[[k]] <- contrast * c(β1, β2) # sparse case
    Prxy$δst_nsprs[[k]] <- contrast # non-sparse case
  }
  b <- ceiling(.12 * p - 1) # p > 8
  Trgt$βt <- rep(c(1, rep(0, p / b)), b)[1:p]
  Trgt$Σξt <- diag(p2)
  return(list(Prxy = Prxy, Trgt = Trgt))
}

# DGP ---------------------------------------------------------------------


dgp.gen.features <- function(n, p1, p2, Σx, Σξ, trgt = F, P = NULL) {
  if (!is.null(Σx)) {
    X <- mvrnorm(n, rep(0, p1), Σx)
  } else {
    X <- matrix(runif(n * p1, -2, 2), n)
    }
    tb <- tibble(X = X)
    if (p2 > 0) {
      if (is.null(P)) {
        ## Non-linear feature mapping inspired by
        ## Additive map in Zhang and Simon, 2023
        fmap.hj <- function(x, j, A) {
          y <- 0
          for (i in A) {
            y <- y + as.numeric(i %% 2 == 1) * (0.5 - abs(x[i] - 0.5)) +
              as.numeric(i %% 2 == 0) * exp(sample(c(-1, 1) / 5, 1) * x[i])
          }
          y
        }
        fmap.h <- function(x, p2, A) {
          yv <- numeric(p2)
          for (j in 1:p2) {
            yv[j] <- fmap.hj(x, j, A[, j])
          }
          yv
        }
        H <- function(X, p2, trgt, A) {
          H <- apply(X, 1, fmap.h, p2 = p2, A = A) %>% t
          if (!trgt) {
            return(H + cos(H))
          }
          H
        }
        A <- replicate(p2, 1:5)
        Z <- H(X, p2, trgt, A)
      } else {
        Z <- X %*% P 
      }
      Ξ <- mvrnorm(n, rep(0, p2), Σξ)
      tb$Z <- Z + Ξ
    }
  return(tb)
}

### Data-Generating Process
# np: sample size of proxy studies (> d, set to be identical across K studies)
# nt: sample size of target study (typically nt{\gg}log(p1 + p2))
# p1, p2: # of matched and mismatched features (may depend on nt)
# sparse: whether beta contrast should be sparse (proximity)
dgp.gen.D = function(seedno, fmap, np, nt, ntest, Prxy, Trgt) {
  set.seed(seedno)
  K <- Prxy$Σ1p %>% length
  p1 <- Prxy$Σ1p[[1]] %>% ncol
  p2 <- Prxy$Σξp[[1]] %>% ncol
  p <- p1 + p2
  𝒫 <- list()
  δst <- list()
  E.cum <- mvrnorm(max(np * K, nt + ntest), c(0, 0), matrix(c(1, .5, .5, 1), 2))
  
  ## Target domain
  Pt <- NULL
  if (fmap == 'linear') Pt <- matrix(rbeta(p1 * p2, 10, 10), p1) * 10
  εt <- E.cum[1:nt, 1]
  εtest <- E.cum[1:ntest + nt, 1]
  Dtrgt <- dgp.gen.features(nt + ntest, p1, p2, NULL, Trgt$Σξt[1:p2, 1:p2], 
                            P = Pt, trgt = T)
  Dt <- Dtrgt[1:nt,]
  Dtest <- Dtrgt[1:ntest + nt,]
  Yt <- as.matrix(Dt) %*% Trgt$βt[1:p] + εt
  Ytest <- as.matrix(Dtest) %*% Trgt$βt[1:p] + εtest
  
  𝒯 <- tibble(Y = c(Yt), X = Dt$X) # missing Z
  testdata = tibble(Y = c(Ytest), X = Dtest$X, εtest = εtest)
  
  ## Proxy domain
  for (k in 1:K) {
    εp <- E.cum[(np * (k - 1) + 1):(np * k), 2]
    Pp <- NULL
    if (fmap == 'linear') Pp <- Pt + 1 / 3 * (matrix(rbeta(p1 * p2, 4, 4), p1) - 1 / 2)
    Dp <- dgp.gen.features(np, p1, p2, Prxy$Σxp[[k]][1:p1, 1:p1], 
                           Prxy$Σξp[[k]][1:p2, 1:p2], P = Pp)
    δst[['sp']] <- Prxy$δst_sprs[[k]][1:p]
    δst[['nsp']] <- Prxy$δst_nsprs[[k]][1:p]
    for (δsti in c('sp', 'nsp')) {
      βpk <- Trgt$βt[1:p] - δst[[δsti]][1:p]
      Yp <- as.matrix(Dp) %*% βpk + εp
      proxydata = tibble(Y = c(Yp), D = Dp)
      𝒫[[δsti]][[k]] <- proxydata
    }
  }
  return(
    list(proxydata = 𝒫, targetdata = 𝒯, testdata = testdata,
      p1 = p1, p2 = p2, K = K, βt = Trgt$βt[1:p])
  )
}