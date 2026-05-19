scale_rows <- function(mat, center_ids = NULL) {
  mat <- as.matrix(mat)
  if (is.null(center_ids)) {
    center_ids <- colnames(mat)
  }
  center_ids <- intersect(center_ids, colnames(mat))
  if (length(center_ids) == 0) {
    stop("No `center_ids` are present in matrix columns.", call. = FALSE)
  }

  centers <- rowMeans(mat[, center_ids, drop = FALSE], na.rm = TRUE)
  scales <- apply(mat[, center_ids, drop = FALSE], 1, stats::sd, na.rm = TRUE)
  scales[is.na(scales) | scales == 0] <- 1
  sweep(sweep(mat, 1, centers, "-"), 1, scales, "/")
}

standardize_sex <- function(x) {
  x <- as.character(x)
  x <- ifelse(x %in% c("m", "M", "male", "Male"), "Male", x)
  x <- ifelse(x %in% c("f", "F", "female", "Female"), "Female", x)
  x
}

fit_loess_predict <- function(x, y, new_x, span = 1.5) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]

  if (length(unique(x)) < 4 || length(y) < 6) {
    fit <- rep(NA_real_, length(new_x))
    se <- rep(NA_real_, length(new_x))
    return(list(fit = fit, se = se))
  }

  data <- data.frame(x = x, y = y)
  model <- stats::loess(
    y ~ x,
    data = data,
    span = span
  )
  pred <- stats::predict(model, newdata = data.frame(x = new_x), se = TRUE)
  list(fit = as.numeric(pred$fit), se = as.numeric(pred$se.fit))
}

resolve_trajectory_span <- function(span, feature, sex, tissue) {
  if (is.numeric(span) && length(span) == 1) {
    return(span)
  }

  if (is.list(span) && !is.data.frame(span)) {
    if (!is.null(span[[tissue]])) {
      return(resolve_trajectory_span(span[[tissue]], feature, sex, tissue))
    }
    stop("Span list must contain entries named `Normal` and/or `Tumor`.", call. = FALSE)
  }

  if (is.data.frame(span)) {
    required <- c("feature", "sex", "tissue", "span")
    missing <- setdiff(required, names(span))
    if (length(missing) > 0) {
      stop(
        "Span data frame is missing required column(s): ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }

    hit <- span$feature == feature & span$sex == sex & span$tissue == tissue
    if (!any(hit)) {
      stop(
        "No span specified for feature=", feature,
        ", sex=", sex, ", tissue=", tissue,
        call. = FALSE
      )
    }
    return(span$span[which(hit)[1]])
  }

  stop(
    "`span` must be a single number, a list, or a data frame with ",
    "feature/sex/tissue/span columns.",
    call. = FALSE
  )
}

adjust_normal_feature <- function(values, metadata, covariates) {
  covariates <- covariates[covariates %in% names(metadata)]
  if (length(covariates) == 0) {
    return(as.numeric(scale(values)))
  }

  data <- data.frame(score = as.numeric(values), metadata[, covariates, drop = FALSE])

  # Match the manuscript trajectory code: impute missing pH to the mean and
  # combine H with C before using C as the ethnicity reference level.
  if ("pH" %in% names(data)) {
    data$pH[is.na(data$pH)] <- mean(data$pH, na.rm = TRUE)
  }
  if ("Ethnicity" %in% names(data)) {
    data$Ethnicity <- gsub("H", "C", as.character(data$Ethnicity))
    data$Ethnicity <- stats::relevel(factor(data$Ethnicity), ref = "C")
  }

  data <- data[stats::complete.cases(data), , drop = FALSE]
  if (nrow(data) < length(covariates) + 3) {
    return(rep(NA_real_, length(values)))
  }

  for (covariate in covariates) {
    if (is.character(data[[covariate]])) {
      data[[covariate]] <- factor(data[[covariate]])
    }
  }

  formula <- stats::as.formula(paste("score ~", paste(covariates, collapse = " + ")))
  fit <- stats::lm(formula, data = data)

  model_metadata <- data.frame(metadata[, covariates, drop = FALSE])
  if ("pH" %in% names(model_metadata)) {
    model_metadata$pH[is.na(model_metadata$pH)] <- mean(model_metadata$pH, na.rm = TRUE)
  }
  if ("Ethnicity" %in% names(model_metadata)) {
    model_metadata$Ethnicity <- gsub("H", "C", as.character(model_metadata$Ethnicity))
    model_metadata$Ethnicity <- factor(model_metadata$Ethnicity, levels = levels(data$Ethnicity))
  }

  out <- rep(NA_real_, length(values))
  usable <- stats::complete.cases(model_metadata) & !is.na(values)
  pred <- stats::predict(fit, newdata = model_metadata[usable, , drop = FALSE])
  out[usable] <- values[usable] - pred
  as.numeric(scale(out))
}

#' Compare normal and tumor age trajectories
#'
#' This function implements the Figure 2-style normal/tumor trajectory
#' comparison using public tumor molecular data and package-stored normal
#' reference data. Normal values are first adjusted for technical covariates
#' such as pH, PMI, and ethnicity, then normal and tumor age trajectories are
#' fit separately by sex.
#'
#' Manuscript reproducibility details are intentionally preserved. The loess
#' fits use base [stats::loess()] defaults, matching the original trajectory
#' scripts. For normal/reference data, missing pH should be imputed before
#' sex-stratified fitting, `Ethnicity == "H"` is combined with `"C"`, and
#' `"C"` is used as the reference level when the covariate model includes
#' ethnicity. The Figure 2C reproduction script also filters tumor fitting
#' samples to `age <= 80`, evaluates each sex-stratified curve at all tumor
#' ages `<= 50`, and uses manuscript-specific adaptive spans.
#'
#' @param tumor_mat Tumor feature-by-sample matrix.
#' @param tumor_metadata Tumor sample metadata.
#' @param normal_mat Normal/reference feature-by-sample matrix.
#' @param normal_metadata Normal/reference sample metadata.
#' @param features Features to plot/model.
#' @param tumor_sample_col Tumor metadata sample ID column.
#' @param tumor_age_col Tumor metadata age column.
#' @param tumor_sex_col Tumor metadata sex column.
#' @param normal_sample_col Normal metadata sample ID column.
#' @param normal_age_col Normal metadata age column.
#' @param normal_sex_col Normal metadata sex column.
#' @param normal_covariates Covariates used to adjust normal data.
#' @param center_age_range Age range used to center/scale tumor and normal data.
#' @param span Loess span for trajectory fitting. This can be a single numeric
#'   value used for all fits, a list with `Normal` and `Tumor` entries, or a
#'   data frame with columns `feature`, `sex`, `tissue`, and `span`.
#' @param prediction_ages Optional numeric vector of ages where trajectories
#'   should be predicted. If `NULL`, trajectories are predicted at the tumor
#'   sample ages for the sex being modeled. For manuscript Figure 2C
#'   reproduction, pass all tumor sample ages `<= 50` so each sex-stratified
#'   curve is evaluated on the same age support used by the original `tn.df`
#'   object.
#' @param ci_level Confidence level for fitted trajectories.
#'
#' @return A data frame with fitted values, standard errors, and confidence
#'   intervals for normal and tumor trajectories at tumor sample ages or at
#'   `prediction_ages` when provided.
#' @export
ageTMP_compare_normal_tumor_trajectory <- function(
  tumor_mat,
  tumor_metadata,
  normal_mat,
  normal_metadata,
  features,
  tumor_sample_col = "id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  normal_sample_col = "ID",
  normal_age_col = "Age",
  normal_sex_col = "Gender",
  normal_covariates = c("pH", "PMI", "Ethnicity"),
  center_age_range = c(0, 26),
  span = 1.5,
  prediction_ages = NULL,
  ci_level = 0.95
) {
  tumor_mat <- as.matrix(tumor_mat)
  normal_mat <- as.matrix(normal_mat)
  features <- intersect(features, intersect(rownames(tumor_mat), rownames(normal_mat)))
  if (length(features) == 0) {
    stop("No requested features are present in both tumor and normal matrices.", call. = FALSE)
  }

  tumor_metadata <- data.frame(tumor_metadata, stringsAsFactors = FALSE)
  normal_metadata <- data.frame(normal_metadata, stringsAsFactors = FALSE)
  tumor_metadata[[tumor_sample_col]] <- ageTMP_normalize_sample_ids(tumor_metadata[[tumor_sample_col]])
  tumor_metadata[[tumor_sex_col]] <- standardize_sex(tumor_metadata[[tumor_sex_col]])
  normal_metadata[[normal_sex_col]] <- standardize_sex(normal_metadata[[normal_sex_col]])

  # Match the manuscript code: pH imputation is done once on the full
  # normal/reference metadata before sex-stratified modeling.
  normal_metadata_model <- normal_metadata
  if ("pH" %in% names(normal_metadata_model)) {
    normal_metadata_model$pH[is.na(normal_metadata_model$pH)] <-
      mean(normal_metadata_model$pH, na.rm = TRUE)
  }

  colnames(tumor_mat) <- ageTMP_normalize_sample_ids(colnames(tumor_mat))
  tumor_ids <- intersect(colnames(tumor_mat), tumor_metadata[[tumor_sample_col]])
  tumor_metadata <- tumor_metadata[match(tumor_ids, tumor_metadata[[tumor_sample_col]]), , drop = FALSE]
  tumor_mat <- tumor_mat[features, tumor_ids, drop = FALSE]

  normal_ids <- intersect(colnames(normal_mat), normal_metadata[[normal_sample_col]])
  normal_metadata_match <- match(normal_ids, normal_metadata[[normal_sample_col]])
  normal_metadata <- normal_metadata[normal_metadata_match, , drop = FALSE]
  normal_metadata_model <- normal_metadata_model[normal_metadata_match, , drop = FALSE]
  normal_mat <- normal_mat[features, normal_ids, drop = FALSE]

  center_ids <- tumor_metadata[[tumor_sample_col]][
    tumor_metadata[[tumor_age_col]] >= center_age_range[1] &
      tumor_metadata[[tumor_age_col]] <= center_age_range[2]
  ]
  tumor_scaled <- scale_rows(tumor_mat, center_ids = center_ids)

  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  out <- list()

  for (feature in features) {
    for (sex in c("Male", "Female")) {
      tumor_sex_idx <- which(tumor_metadata[[tumor_sex_col]] == sex)
      normal_sex_idx <- which(normal_metadata[[normal_sex_col]] == sex)
      if (length(tumor_sex_idx) == 0 || length(normal_sex_idx) == 0) {
        next
      }

      tumor_age <- tumor_metadata[[tumor_age_col]][tumor_sex_idx]
      tumor_values <- as.numeric(tumor_scaled[feature, tumor_sex_idx])
      tumor_span <- resolve_trajectory_span(span, feature, sex, "Tumor")
      normal_span <- resolve_trajectory_span(span, feature, sex, "Normal")

      if (is.null(prediction_ages)) {
        pred_age <- tumor_age
        pred_sample_id <- tumor_metadata[[tumor_sample_col]][tumor_sex_idx]
        pred_observed <- tumor_values
      } else {
        pred_age <- sort(unique(as.numeric(prediction_ages)))
        pred_sample_id <- rep(NA_character_, length(pred_age))
        pred_observed <- rep(NA_real_, length(pred_age))
      }

      tumor_pred <- fit_loess_predict(tumor_age, tumor_values, pred_age, span = tumor_span)

      normal_values <- adjust_normal_feature(
        as.numeric(normal_mat[feature, normal_sex_idx]),
        normal_metadata_model[normal_sex_idx, , drop = FALSE],
        normal_covariates
      )
      normal_age <- normal_metadata[[normal_age_col]][normal_sex_idx]
      normal_center <- normal_age >= center_age_range[1] & normal_age <= center_age_range[2]
      normal_values <- normal_values - mean(normal_values[normal_center], na.rm = TRUE)
      normal_sd <- stats::sd(normal_values[normal_center], na.rm = TRUE)
      if (!is.na(normal_sd) && normal_sd > 0) {
        normal_values <- normal_values / normal_sd
      }
      normal_pred <- fit_loess_predict(normal_age, normal_values, pred_age, span = normal_span)

      base <- data.frame(
        feature = feature,
        sample_id = pred_sample_id,
        age = pred_age,
        sex = sex,
        stringsAsFactors = FALSE
      )

      out[[length(out) + 1]] <- data.frame(
        base,
        tissue = "Tumor",
        observed = pred_observed,
        fit = tumor_pred$fit,
        se = tumor_pred$se,
        ci_lower = tumor_pred$fit - z * tumor_pred$se,
        ci_upper = tumor_pred$fit + z * tumor_pred$se,
        span = tumor_span,
        stringsAsFactors = FALSE
      )
      out[[length(out) + 1]] <- data.frame(
        base,
        tissue = "Normal",
        observed = NA_real_,
        fit = normal_pred$fit,
        se = normal_pred$se,
        ci_lower = normal_pred$fit - z * normal_pred$se,
        ci_upper = normal_pred$fit + z * normal_pred$se,
        span = normal_span,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, out)
}
