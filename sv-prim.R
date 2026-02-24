## additive map by Zhang and Simon, 2023
hj <- function(x, j, A) {
  y <- 0
  for (i in A) {
    y <- y + as.numeric(i %% 2 == 1) * (0.5 - abs(x[i] - 0.5)) +
      as.numeric(i %% 2 == 0) * exp(-x[i])
  }
  y
}

h <- function(x, p2, A) {
  yv <- numeric(p2)
  for (j in 1:p2) {
    yv[j] <- hj(x, j, A[, j])
  }
  yv
}

H <- function(X, p2, A, trgt) {
  H <- apply(X, 1, h, p2 = p2, A = A) %>% t
  if (!trgt) {
    return(H + sin(H))
  }
  H
}