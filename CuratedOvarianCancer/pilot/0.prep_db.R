# Load curated ovarian cancer datasets into analysis-ready matrices.
# Follows preprocessing in ../1.data_pre_processing.R:
#   - expression: samples x genes (transpose of ExpressionSet assay matrix)
#   - tumor samples when sample_type is annotated
#   - analysis subset: vital_status == "deceased" with valid days_to_death
#   - outcome: y = log(days_to_death)
#
# Dataset selection rules are documented in DATASET_SELECTION.md.
library(curatedOvarianData)
library(limma)

PILOT_SELECTION_CRITERIA <- list(
  min_number_of_events = 15L,
  min_sample_size = 40L,
  min_number_of_genes = 1000L,
  tumor_only = TRUE
)

DATASET_EXCLUSION_REGISTRY <- data.frame(
  id = c(
    "GSE12418", "GSE12470", "GSE20565", "GSE2109", "GSE44104",
    "GSE6008", "GSE6822", "PMID15897565",
    "GSE30009", "TCGA.mirna.8x15kv2",
    "GSE32063", "TCGA.RNASeqV2",
    "GSE8842", "PMID19318476", "GSE19829.GPL570"
  ),
  stage = c(
    rep("hard", 8L),
    rep("modality", 2L),
    rep("duplicate_cohort", 2L),
    rep("quality", 3L)
  ),
  reason = c(
    "No vital_status / days_to_death survival annotation",
    "No survival annotation; includes healthy controls",
    "Differential diagnosis study; no survival endpoint",
    "IGC multi-tumor atlas; no survival endpoint",
    "No survival annotation",
    "In vitro LPA-stimulation profile; no survival endpoint",
    "Tumor classification study; no survival endpoint; includes cell lines",
    "No survival annotation in curated metadata",
    "Pre-defined 359-gene signature panel (p < 1,000)",
    "miRNA platform; incompatible with mRNA transfer learning",
    "Duplicate of GSE32062.GPL6480 (smaller validation subset)",
    "Duplicate TCGA cohort on RNA-seq platform; keep TCGA microarray",
    "Fewer than 15 deceased events among tumor samples",
    "Early-stage cohort; excluded by vignette createEsetList.R",
    "Fewer than 40 tumor samples"
  ),
  stringsAsFactors = FALSE
)
#' List all ExpressionSet objects shipped with curatedOvarianData.
list_curated_ovarian_datasets <- function() {
  pkg_info <- utils::data(
    package = "curatedOvarianData",
    envir = environment(),
    verbose = FALSE
  )
  eset_names <- sub(" .*", "", pkg_info$results[, "Item"])
  sort(unique(eset_names))
}
eset_name_from_id <- function(dataset_id) {
  if (grepl("_eset$", dataset_id)) {
    return(dataset_id)
  }
  paste0(dataset_id, "_eset")
}
dataset_id_from_eset <- function(eset_name) {
  sub("_eset$", "", eset_name)
}
#' Return exclusion table for datasets not in the pilot pool.
get_dataset_exclusions <- function() {
  DATASET_EXCLUSION_REGISTRY
}
#' Evaluate whether a loaded dataset object passes pilot quality filters.
passes_pilot_criteria <- function(
    ds,
    criteria = PILOT_SELECTION_CRITERIA) {
  if (!is.null(ds$load_error)) {
    return(FALSE)
  }
  ds$n_genes >= criteria$min_number_of_genes &&
    ds$n_tumor >= criteria$min_sample_size &&
    ds$n_analysis >= criteria$min_number_of_events
}
#' Filter eset names to the pilot inclusion pool.
filter_pilot_eset_names <- function(
    eset_names = list_curated_ovarian_datasets(),
    criteria = PILOT_SELECTION_CRITERIA,
    exclusion_registry = DATASET_EXCLUSION_REGISTRY) {
  excluded_ids <- exclusion_registry$id
  candidate_names <- eset_names[
    !dataset_id_from_eset(eset_names) %in% excluded_ids
  ]
  eligible <- vapply(candidate_names, function(eset_name) {
    ds <- eset_to_dataset(
      eset_name,
      deceased_only = TRUE,
      tumor_only = criteria$tumor_only
    )
    passes_pilot_criteria(ds, criteria = criteria)
  }, logical(1))
  sort(candidate_names[eligible])
}
#' Convert one ExpressionSet to an analysis list element.
#'
#' @param eset_name Character name of the data object (e.g. "GSE9891_eset").
#' @param deceased_only If TRUE, restrict X and y to deceased subjects with valid outcome.
#' @param tumor_only If TRUE, keep tumor samples when sample_type is annotated.
eset_to_dataset <- function(
    eset_name,
    deceased_only = TRUE,
    tumor_only = TRUE) {
  data(
    list = eset_name,
    package = "curatedOvarianData",
    envir = environment()
  )
  eset <- get(eset_name, envir = environment())
  expr_mat <- t(limma::as.matrix.ExpressionSet(eset))
  pheno <- Biobase::pData(eset)
  sample_ids <- rownames(expr_mat)
  if (is.null(sample_ids)) {
    sample_ids <- colnames(eset)
    rownames(expr_mat) <- sample_ids
  }
  tumor_idx <- seq_len(nrow(expr_mat))
  if (tumor_only && "sample_type" %in% colnames(pheno)) {
    tumor_idx <- which(pheno$sample_type == "tumor")
    expr_mat <- expr_mat[tumor_idx, , drop = FALSE]
    pheno <- pheno[tumor_idx, , drop = FALSE]
    sample_ids <- sample_ids[tumor_idx]
  }
  dataset_id <- dataset_id_from_eset(eset_name)
  pkg_row <- utils::data(
    package = "curatedOvarianData",
    envir = environment(),
    verbose = FALSE
  )$results
  title_row <- pkg_row[pkg_row[, "Item"] == eset_name, "Title", drop = TRUE]
  has_vital <- "vital_status" %in% colnames(pheno)
  has_days <- "days_to_death" %in% colnames(pheno)
  deceased_idx <- if (has_vital) {
    which(pheno$vital_status == "deceased")
  } else {
    integer()
  }
  days_to_death <- if (has_days) {
    as.numeric(pheno$days_to_death)
  } else {
    rep(NA_real_, nrow(expr_mat))
  }
  valid_outcome_idx <- which(!is.na(days_to_death) & days_to_death > 0)
  analysis_idx <- intersect(deceased_idx, valid_outcome_idx)
  y_all <- rep(NA_real_, nrow(expr_mat))
  y_all[valid_outcome_idx] <- log(days_to_death[valid_outcome_idx])
  out <- list(
    id = dataset_id,
    eset_name = eset_name,
    title = if (length(title_row)) unname(title_row) else NA_character_,
    expression = expr_mat,
    pheno = pheno,
    sample_ids = sample_ids,
    genes = colnames(expr_mat),
    tumor_idx = seq_len(nrow(expr_mat)),
    deceased_idx = deceased_idx,
    valid_outcome_idx = valid_outcome_idx,
    analysis_idx = analysis_idx,
    y_all = y_all,
    n_samples = nrow(expr_mat),
    n_tumor = nrow(expr_mat),
    n_genes = ncol(expr_mat),
    n_deceased = length(deceased_idx),
    n_valid_outcome = length(valid_outcome_idx),
    n_analysis = length(analysis_idx),
    has_vital_status = has_vital,
    has_days_to_death = has_days,
    tumor_only = tumor_only,
    analysis_ready = length(analysis_idx) >= 2L
  )
  if (deceased_only) {
    out$X <- expr_mat[analysis_idx, , drop = FALSE]
    out$y <- y_all[analysis_idx]
    out$sample_ids_analysis <- sample_ids[analysis_idx]
  } else {
    out$X <- expr_mat
    out$y <- y_all
    out$sample_ids_analysis <- sample_ids
  }
  out
}
#' Load every available dataset and return a named list of analysis objects.
load_all_ovarian_datasets <- function(
    eset_names = NULL,
    deceased_only = TRUE,
    tumor_only = TRUE,
    verbose = TRUE) {
  if (is.null(eset_names)) {
    eset_names <- list_curated_ovarian_datasets()
  }
  datasets <- setNames(
    vector("list", length(eset_names)),
    dataset_id_from_eset(eset_names)
  )
  for (eset_name in eset_names) {
    ds <- tryCatch(
      eset_to_dataset(
        eset_name,
        deceased_only = deceased_only,
        tumor_only = tumor_only
      ),
      error = function(e) {
        list(
          id = dataset_id_from_eset(eset_name),
          eset_name = eset_name,
          load_error = conditionMessage(e),
          analysis_ready = FALSE
        )
      }
    )
    datasets[[ds$id]] <- ds
    if (verbose) {
      if (!is.null(ds$load_error)) {
        message(ds$id, ": FAILED (", ds$load_error, ")")
      } else {
        message(
          ds$id,
          ": n_tumor=", ds$n_tumor,
          ", p=", ds$n_genes,
          ", deceased=", ds$n_deceased,
          ", analysis=", ds$n_analysis,
          if (!ds$analysis_ready) " [not analysis-ready]" else ""
        )
      }
    }
  }
  attr(datasets, "n_total") <- length(datasets)
  attr(datasets, "n_analysis_ready") <- sum(vapply(
    datasets,
    function(ds) isTRUE(ds$analysis_ready),
    logical(1)
  ))
  datasets
}
#' Load datasets selected for the pilot study.
load_pilot_ovarian_datasets <- function(
    criteria = PILOT_SELECTION_CRITERIA,
    verbose = TRUE) {
  eset_names <- filter_pilot_eset_names(criteria = criteria)
  datasets <- load_all_ovarian_datasets(
    eset_names = eset_names,
    deceased_only = TRUE,
    tumor_only = criteria$tumor_only,
    verbose = verbose
  )
  attr(datasets, "selection_criteria") <- criteria
  attr(datasets, "target_id") <- "GSE9891"
  attr(datasets, "reference_proxy_id") <- "GSE26712"
  datasets
}
#' Reduce one dataset object to regression-ready matrices.
#'
#' @return A list with \code{X} (n x p) and \code{y} (length n).
as_regression_data <- function(ds) {
  if (is.null(ds$X) || is.null(ds$y)) {
    stop("Dataset ", ds$id %||% "<unknown>", " is missing X or y.", call. = FALSE)
  }
  list(
    X = as.matrix(ds$X),
    y = as.numeric(ds$y)
  )
}
#' Convert a named dataset list to regression-ready form.
#'
#' Each element is a list \code{list(X = ..., y = ...)} suitable for direct use
#' in \code{glmnet(X, y)} and the pilot HTL / SAND pipeline.
datasets_to_regression_list <- function(datasets) {
  ids <- names(datasets)
  reg_list <- setNames(
    lapply(datasets, as_regression_data),
    ids
  )
  for (id in ids) {
    n <- nrow(reg_list[[id]]$X)
    if (length(reg_list[[id]]$y) != n) {
      stop(
        "Length mismatch for ", id, ": nrow(X) = ", n,
        ", length(y) = ", length(reg_list[[id]]$y),
        call. = FALSE
      )
    }
  }
  attrs <- c(
    "selection_criteria", "target_id", "reference_proxy_id",
    "n_total", "n_analysis_ready"
  )
  for (nm in attrs) {
    if (!is.null(attr(datasets, nm))) {
      attr(reg_list, nm) <- attr(datasets, nm)
    }
  }
  reg_list
}
#' Load pilot datasets as a regression-ready named list.
load_pilot_regression_data <- function(
    criteria = PILOT_SELECTION_CRITERIA,
    verbose = TRUE) {
  datasets_to_regression_list(
    load_pilot_ovarian_datasets(criteria = criteria, verbose = verbose)
  )
}
#' Tabular summary of loaded datasets.
summarize_datasets <- function(datasets) {
  ids <- names(datasets)
  do.call(
    rbind,
    lapply(ids, function(id) {
      ds <- datasets[[id]]
      if (!is.null(ds$X) && !is.null(ds$y)) {
        return(data.frame(
          id = id,
          eset_name = ds$eset_name %||% paste0(id, "_eset"),
          n = nrow(ds$X),
          p = ncol(ds$X),
          n_tumor = ds$n_tumor %||% NA_integer_,
          n_genes = ds$n_genes %||% ncol(ds$X),
          n_deceased = ds$n_deceased %||% NA_integer_,
          n_analysis = length(ds$y),
          analysis_ready = nrow(ds$X) >= 2L,
          load_error = ds$load_error %||% NA_character_,
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        id = ds$id %||% id,
        eset_name = ds$eset_name %||% NA_character_,
        n = NA_integer_,
        p = NA_integer_,
        n_tumor = ds$n_tumor %||% NA_integer_,
        n_genes = ds$n_genes %||% NA_integer_,
        n_deceased = ds$n_deceased %||% NA_integer_,
        n_analysis = ds$n_analysis %||% NA_integer_,
        analysis_ready = isTRUE(ds$analysis_ready),
        load_error = ds$load_error %||% NA_character_,
        stringsAsFactors = FALSE
      )
    })
  )
}
`%||%` <- function(x, y) if (is.null(x)) y else x
# --- run when executed as a script -----------------------------------------
if (sys.nframe() == 0L) {
  pilot_dir <- file.path("CuratedOvarianCancer", "pilot")
  datasets_full <- load_pilot_ovarian_datasets()
  ovarian_datasets <- datasets_to_regression_list(datasets_full)
  dataset_summary <- summarize_datasets(ovarian_datasets)
  out_rds <- file.path(pilot_dir, "ovarian_datasets.rds")
  saveRDS(ovarian_datasets, out_rds)
  out_csv <- file.path(pilot_dir, "dataset_summary.csv")
  utils::write.csv(dataset_summary, out_csv, row.names = FALSE)
  out_excl <- file.path(pilot_dir, "dataset_exclusions.csv")
  utils::write.csv(get_dataset_exclusions(), out_excl, row.names = FALSE)
  cat(
    "\nLoaded", attr(ovarian_datasets, "n_total"), "pilot datasets;",
    attr(ovarian_datasets, "n_analysis_ready"), "analysis-ready.\n"
  )
  cat("Target:", attr(ovarian_datasets, "target_id"), "\n")
  cat("Reference proxy:", attr(ovarian_datasets, "reference_proxy_id"), "\n")
  cat("Saved:", out_rds, "\n")
  cat("Saved:", out_csv, "\n")
  cat("Saved:", out_excl, "\n")
}