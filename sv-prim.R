# Add two matrices with zero-padding (no name handling)
add_mats_pad <- function(A, B) {
  if (is.null(A)) return(B)
  R <- max(nrow(A), nrow(B))
  C <- max(ncol(A), ncol(B))
  pad <- function(M) {
    out <- matrix(0, R, C)
    out[seq_len(nrow(M)), seq_len(ncol(M))] <- M
    out
  }
  pad(A) + pad(B)
}

# For cbind with padding
cbind_pad0 <- function(..., fill = 0) {
  xs <- list(...)
  # make sure vectors become column matrices
  xs <- lapply(xs, function(x) {
    if (is.null(dim(x))) matrix(x, ncol = 1) else x
  })
  maxr <- max(sapply(xs, nrow))
  xs_pad <- lapply(xs, function(y) {
    if (nrow(y) < maxr) {
      rbind(y, matrix(fill, nrow = maxr - nrow(y), ncol = ncol(y),
                      dimnames = list(NULL, colnames(y))))
    } else y
  })
  do.call(cbind, xs_pad)
}

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