select_adaptive_span <- function(x, y, min_span = 0.5, max_span = 3, span_step = 0.1) {
  if (!requireNamespace("locfit", quietly = TRUE)) {
    stop("Package `locfit` is required when `adaptive_span = TRUE`.", call. = FALSE)
  }
  suppressPackageStartupMessages(
    require("locfit", quietly = TRUE, character.only = TRUE)
  )

  ok <- !is.na(x) & !is.na(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  if (length(unique(x)) < 4 || length(y) < 6) {
    return(NA_real_)
  }

  spans <- seq(min_span, max_span, by = span_step)
  data <- data.frame(x = x, y = y)
  scores <- vapply(spans, function(alpha) {
    score <- try(locfit::gcv(y ~ x, data = data, alpha = alpha), silent = TRUE)
    if (inherits(score, "try-error")) {
      return(Inf)
    }
    round(as.numeric(score["gcv"]), 3)
  }, numeric(1))

  if (all(!is.finite(scores))) {
    return(NA_real_)
  }
  spans[which.min(scores)]
}

#' Fit tumor age-dependent molecular trajectories
#'
#' Fit sex-stratified age-dependent tumor molecular trajectories from a public
#' feature-by-sample matrix and sample metadata. This is the tumor-only AD-TMP
#' fitting step used by manuscript heatmap-style trajectory panels such as
#' Figure 2A.
#'
#' @details
#' The function intentionally mirrors key details from the manuscript tumor
#' trajectory workflow. Tumor sample IDs are harmonized with
#' [ageTMP_normalize_sample_ids()], rows are centered/scaled using the requested
#' `center_age_range`, optional fitting can be restricted with `fit_age_range`,
#' and sex-stratified [stats::loess()] models are fit with base loess defaults.
#'
#' For exact manuscript reproduction, the original protein workflow used
#' feature/sex-specific adaptive loess spans selected by generalized
#' cross-validation over a span grid. Set `adaptive_span = TRUE` to recompute
#' those spans from the public tumor matrix; this requires the suggested
#' `locfit` package and can be slow for thousands of proteins. A single numeric
#' `span`, or a data frame with `feature`, `sex`, `tissue`, and `span` columns,
#' can also be supplied when spans are known or when a faster reconstruction is
#' desired.
#'
#' @param tumor_mat Tumor feature-by-sample matrix.
#' @param tumor_metadata Tumor sample metadata.
#' @param features Features to model. Defaults to all matrix rows.
#' @param tumor_sample_col Tumor metadata sample ID column.
#' @param tumor_age_col Tumor metadata age column.
#' @param tumor_sex_col Tumor metadata sex column.
#' @param center_age_range Age range used to center and scale each feature.
#' @param fit_age_range Optional age range used for loess fitting.
#' @param span Loess span. May be a single number or a span data frame accepted
#'   by the normal/tumor trajectory functions.
#' @param adaptive_span Recompute feature/sex-specific spans by GCV.
#' @param min_span Minimum span for adaptive selection.
#' @param max_span Maximum span for adaptive selection.
#' @param span_step Span grid step for adaptive selection.
#' @param prediction_ages Optional ages where trajectories are predicted. If
#'   `NULL`, predictions are made at all harmonized tumor sample ages.
#' @param prediction_sample_ids Optional IDs corresponding to `prediction_ages`.
#' @param ci_level Confidence level for fitted trajectories.
#'
#' @return A data frame with one row per feature, sex, and prediction age.
#' @export
ageTMP_fit_tumor_trajectory <- function(
  tumor_mat,
  tumor_metadata,
  features = rownames(tumor_mat),
  tumor_sample_col = "id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  center_age_range = c(0, 26),
  fit_age_range = NULL,
  span = 1.5,
  adaptive_span = FALSE,
  min_span = 0.5,
  max_span = 3,
  span_step = 0.1,
  prediction_ages = NULL,
  prediction_sample_ids = NULL,
  ci_level = 0.95
) {
  tumor_mat <- as.matrix(tumor_mat)
  features <- intersect(features, rownames(tumor_mat))
  if (length(features) == 0) {
    stop("No requested features are present in `tumor_mat`.", call. = FALSE)
  }

  required_cols <- c(tumor_sample_col, tumor_age_col, tumor_sex_col)
  missing_cols <- setdiff(required_cols, names(tumor_metadata))
  if (length(missing_cols) > 0) {
    stop(
      "Tumor metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  tumor_metadata <- data.frame(tumor_metadata, stringsAsFactors = FALSE)
  tumor_metadata[[tumor_sample_col]] <- ageTMP_normalize_sample_ids(tumor_metadata[[tumor_sample_col]])
  tumor_metadata[[tumor_sex_col]] <- standardize_sex(tumor_metadata[[tumor_sex_col]])
  tumor_metadata[[tumor_age_col]] <- as.numeric(tumor_metadata[[tumor_age_col]])
  colnames(tumor_mat) <- ageTMP_normalize_sample_ids(colnames(tumor_mat))

  tumor_ids <- intersect(colnames(tumor_mat), tumor_metadata[[tumor_sample_col]])
  if (length(tumor_ids) == 0) {
    stop("No tumor matrix columns match tumor metadata sample IDs.", call. = FALSE)
  }

  tumor_metadata <- tumor_metadata[match(tumor_ids, tumor_metadata[[tumor_sample_col]]), , drop = FALSE]
  tumor_mat <- tumor_mat[features, tumor_ids, drop = FALSE]

  fit_keep <- !is.na(tumor_metadata[[tumor_age_col]]) & tumor_metadata[[tumor_sex_col]] %in% c("Male", "Female")
  if (!is.null(fit_age_range)) {
    fit_keep <- fit_keep &
      tumor_metadata[[tumor_age_col]] >= fit_age_range[1] &
      tumor_metadata[[tumor_age_col]] <= fit_age_range[2]
  }
  if (!any(fit_keep)) {
    stop("No tumor samples remain after age/sex filtering.", call. = FALSE)
  }

  fit_metadata <- tumor_metadata[fit_keep, , drop = FALSE]
  fit_mat <- tumor_mat[, fit_keep, drop = FALSE]

  center_keep <- fit_metadata[[tumor_age_col]] >= center_age_range[1] &
    fit_metadata[[tumor_age_col]] <= center_age_range[2]
  center_ids <- fit_metadata[[tumor_sample_col]][center_keep]
  tumor_scaled <- scale_rows(fit_mat, center_ids = center_ids)

  if (is.null(prediction_ages)) {
    pred_age <- tumor_metadata[[tumor_age_col]]
    pred_sample_id <- tumor_metadata[[tumor_sample_col]]
  } else {
    pred_age <- as.numeric(prediction_ages)
    if (is.null(prediction_sample_ids)) {
      pred_sample_id <- rep(NA_character_, length(pred_age))
    } else {
      if (length(prediction_sample_ids) != length(pred_age)) {
        stop("`prediction_sample_ids` must have the same length as `prediction_ages`.", call. = FALSE)
      }
      pred_sample_id <- ageTMP_normalize_sample_ids(prediction_sample_ids)
    }
  }

  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  out <- vector("list", length(features) * 2)
  out_i <- 0L

  for (feature in features) {
    for (sex in c("Male", "Female")) {
      sex_idx <- which(fit_metadata[[tumor_sex_col]] == sex)
      if (length(sex_idx) == 0) {
        next
      }
      tumor_age <- fit_metadata[[tumor_age_col]][sex_idx]
      tumor_values <- as.numeric(tumor_scaled[feature, sex_idx])

      tumor_span <- if (isTRUE(adaptive_span)) {
        select_adaptive_span(tumor_age, tumor_values, min_span, max_span, span_step)
      } else {
        resolve_trajectory_span(span, feature, sex, "Tumor")
      }

      pred <- if (is.na(tumor_span)) {
        list(fit = rep(NA_real_, length(pred_age)), se = rep(NA_real_, length(pred_age)))
      } else {
        fit_loess_predict(tumor_age, tumor_values, pred_age, span = tumor_span)
      }

      out_i <- out_i + 1L
      out[[out_i]] <- data.frame(
        feature = feature,
        sample_id = pred_sample_id,
        age = pred_age,
        sex = sex,
        tissue = "Tumor",
        fit = pred$fit,
        se = pred$se,
        ci_lower = pred$fit - z * pred$se,
        ci_upper = pred$fit + z * pred$se,
        span = tumor_span,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, out[seq_len(out_i)])
}
