# Primitive Functions for non-linear feature map --------------------------

.calc_rmse <- function(y_trth, y_prd) {
  if(length(y_trth) != length(y_prd)) print(mget(ls()))
  sqrt(mean((y_trth - y_prd)^2))
}

.index_spliter <- function(array, n_folds = 5) {
  len <- length(array)
  fold_id = sample(rep(1:n_folds, length.out = len))
  split(seq_len(len), fold_id)
}

# cbinds vectors of varying lengths into a matrix, padding with zeros as needed
# each vector is of class dgCMatrix
cbind_padded_vecs <- function(vec_list) {
  max_len <- max(sapply(vec_list, length))
  padded_mats <- lapply(vec_list, function(vec) {
    if (length(vec) < max_len) {
      padding <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), 
                              dims = c(max_len - length(vec), 1))
      return(rbind(vec, padding))
    } else {
      return(vec)
    }
  })
  do.call(cbind, padded_mats)
}

sum_padded_sparse_matrices <- function(
  mat_list,
  mode = c("auto", "sparse", "dense"),
  density_threshold = 0.10,
  size_threshold = 5e5,
  return_sparse = TRUE
) {
  mode <- match.arg(mode)
  if (length(mat_list) == 0) stop("empty list.")
  max_r <- max(vapply(mat_list, nrow, numeric(1)))
  max_c <- max(vapply(mat_list, ncol, numeric(1)))
  total_cells <- sum(vapply(mat_list, function(mat) nrow(mat) * ncol(mat), numeric(1)))
  total_nnz <- sum(vapply(mat_list, function(mat) {
    if (inherits(mat, "sparseMatrix")) {
      nnzero(mat)
    } else {
      sum(mat != 0)
    }
  }, numeric(1)))
  observed_density <- if (total_cells == 0) 0 else total_nnz / total_cells

  # For dense-ish and small enough problems, direct dense accumulation is faster.
  use_dense <- switch(
    mode,
    dense = TRUE,
    sparse = FALSE,
    auto = (observed_density >= density_threshold) && ((max_r * max_c) <= size_threshold)
  )

  if (use_dense) {
    result_dense <- matrix(0, nrow = max_r, ncol = max_c)
    for (k in seq_along(mat_list)) {
      mat_k <- mat_list[[k]]
      nr <- nrow(mat_k)
      nc <- ncol(mat_k)
      result_dense[seq_len(nr), seq_len(nc)] <- result_dense[seq_len(nr), seq_len(nc)] + as.matrix(mat_k)
    }
    if (isTRUE(return_sparse)) {
      return(Matrix(result_dense, sparse = TRUE))
    }
    return(result_dense)
  }

  n <- length(mat_list)
  i_list <- vector("list", n)
  j_list <- vector("list", n)
  x_list <- vector("list", n)
  for (k in seq_len(n)) {
    # Force general sparse form to preserve all entries (avoid symmetric coercion).
    mat_T <- as(as(as(mat_list[[k]], "dMatrix"), "generalMatrix"), "TsparseMatrix")
    i_list[[k]] <- mat_T@i + 1
    j_list[[k]] <- mat_T@j + 1
    x_list[[k]] <- mat_T@x
  }
  result_mat <- sparseMatrix(
    i = unlist(i_list),
    j = unlist(j_list),
    x = unlist(x_list),
    dims = c(max_r, max_c) # auto Zero-padding
  )
  if (!isTRUE(return_sparse)) {
    return(as.matrix(result_mat))
  }
  return(result_mat)
}