#' Load normal/reference molecular data used for AD-TMP trajectory analysis
#'
#' The normal/reference data used by the paper's age-dependent tumor molecular
#' phenotype (AD-TMP) analyses are distributed with `ageTMP` as a documented RDS
#' object. These data were originally represented in the paper-analysis code by
#' `for_tadj.RData`.
#'
#' The normal developmental reference data come from the DEveLopmental
#' Trajectory Atlas (DELTA) in dorsolateral prefrontal cortex (DLPFC), PMID:
#' 30518843. DELTA is available at <http://amp.pharm.mssm.edu/DELTA>.
#'
#' In the manuscript protein trajectory workflow, these normal/reference data
#' are adjusted with `score ~ pH + PMI + Ethnicity` after imputing missing pH
#' from the full normal-reference metadata, combining ethnicity label `H` with
#' `C`, and using `C` as the reference level. These details are implemented in
#' [ageTMP_compare_normal_tumor_trajectory()] for reproducibility.
#'
#' @param path Optional path to a normal-reference RDS file. If `NULL`, the
#'   package copy at `inst/extdata/normal_reference.rds` is used.
#'
#' @return A named list containing normal/reference matrices, sample metadata,
#'   and provenance metadata.
#' @export
ageTMP_load_normal_reference <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("extdata", "normal_reference.rds", package = "ageTMP")
  }

  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "Normal-reference data were not found. Expected package data at ",
      "`inst/extdata/normal_reference.rds`, or provide `path =` explicitly. ",
      "Build this object from the study normal data previously represented by ",
      "`for_tadj.RData`.",
      call. = FALSE
    )
  }

  readRDS(path)
}
