#' List public ageTMP manuscript data sources
#'
#' @param data_dir Path to the public data directory.
#'
#' @return A data frame with expected source files and whether each exists.
#' @export
ageTMP_data_sources <- function(data_dir = "data") {
  sources <- data.frame(
    source = c(
      "clinical",
      "protein",
      "rna",
      "glyco",
      "phospho",
      "mutation",
      "full_mutation",
      "STable1",
      "STable4",
      "STable5"
    ),
    file = c(
      "clinical_data_04032026.tsv",
      "cDisc_proteome_imputed_data_09152023.tsv",
      "cDisc_rna_coding_10192023.tsv",
      "Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv",
      "cDisc_phosphosite_imputed_data_ischemia_removed_motif_11032023.tsv",
      "cDisc_mutation_10192023.tsv",
      "Disc_full_mutation_data_100224.tsv",
      "STable1.xlsx",
      "STable4.xlsx",
      "STable5.xlsx"
    ),
    stringsAsFactors = FALSE
  )

  sources$path <- file.path(data_dir, sources$file)
  sources$exists <- file.exists(sources$path)
  sources
}

#' Read a public supplementary workbook sheet
#'
#' @param data_dir Path to the public data directory.
#' @param table Supplementary table name, such as `"STable4"`.
#' @param sheet Sheet name.
#' @param ... Additional arguments passed to [readxl::read_excel()].
#'
#' @return A tibble containing the requested sheet.
#' @export
ageTMP_load_supplement <- function(data_dir = "data", table, sheet, ...) {
  if (missing(table) || missing(sheet)) {
    stop("Both `table` and `sheet` are required.", call. = FALSE)
  }

  file <- file.path(data_dir, paste0(table, ".xlsx"))
  if (!file.exists(file)) {
    stop("Supplementary workbook not found: ", file, call. = FALSE)
  }

  sheets <- readxl::excel_sheets(file)
  if (!sheet %in% sheets) {
    stop(
      "Sheet `", sheet, "` not found in ", basename(file),
      ". Available sheets: ", paste(sheets, collapse = ", "),
      call. = FALSE
    )
  }

  readxl::read_excel(file, sheet = sheet, ...)
}

#' Load public cDisc clinical data
#'
#' @param data_dir Path to the public data directory.
#'
#' @return A data frame with clinical metadata.
#' @export
ageTMP_load_clinical <- function(data_dir = "data") {
  file <- file.path(data_dir, "clinical_data_04032026.tsv")
  if (!file.exists(file)) {
    stop("Clinical file not found: ", file, call. = FALSE)
  }

  utils::read.delim(file, sep = "\t", check.names = FALSE)
}

#' Load a public molecular data table
#'
#' @param data_dir Path to the public data directory.
#' @param modality One of `"protein"`, `"rna"`, `"glyco"`, `"phospho"`,
#'   `"mutation"`, or `"full_mutation"`.
#'
#' @return A data frame containing feature annotation columns followed by sample columns.
#' @export
ageTMP_load_molecular <- function(
  data_dir = "data",
  modality = c("protein", "rna", "glyco", "phospho", "mutation", "full_mutation")
) {
  modality <- match.arg(modality)
  files <- c(
    protein = "cDisc_proteome_imputed_data_09152023.tsv",
    rna = "cDisc_rna_coding_10192023.tsv",
    glyco = "Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv",
    phospho = "cDisc_phosphosite_imputed_data_ischemia_removed_motif_11032023.tsv",
    mutation = "cDisc_mutation_10192023.tsv",
    full_mutation = "Disc_full_mutation_data_100224.tsv"
  )

  file <- file.path(data_dir, files[[modality]])
  if (!file.exists(file)) {
    stop("Molecular file not found: ", file, call. = FALSE)
  }

  utils::read.delim(file, sep = "\t", check.names = FALSE)
}
