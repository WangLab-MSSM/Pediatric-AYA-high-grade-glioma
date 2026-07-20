#!/usr/bin/env Rscript

# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai

# -----------------------------------------------------------------------------
# Figure 6N: ITGAV sialylated glycopeptide vs PD-L1/CD274 protein abundance
#
# Nicole L. Tignor, PhD
# Department of Genetics and Genomics
# Icahn School of Medicine at Mount Sinai
#
# Purpose:
#   Reproduce the Figure 6N scatter/correlation panel in a standalone,
#   GitHub-ready form. The original plotting block used hidden in-memory
#   objects (`datac`, `plotme0`, `glyco.data.v2`, `glyco.anno.v2`, `sex.col`).
#   This script replaces those objects with explicit reads from repository
#   tables and keeps the original ggplot/stat_cor/linear-smooth visual logic.
#
# Output:
#   Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation.pdf
# -----------------------------------------------------------------------------

required_packages <- c("ggplot2", "ggpubr", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

library(ggplot2)
library(ggpubr)
library(readxl)




input_stable6 <- file.path("../data", "STable6.xlsx")
input_clinical <- file.path("../data", "STable1.xlsx")
clinical_sheet <- "ClinicalTable"
input_proteome <- file.path("../data", "cDisc_proteome_imputed_data_09152023.tsv")
input_glyco <- file.path("../data", "Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv")

resolve_script_dir <- function(script_name) {
  for (frame in rev(sys.frames())) {
    ofile <- frame$ofile
    if (!is.null(ofile) && basename(ofile) == script_name && file.exists(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)))
  }
  if (file.exists(file.path(getwd(), script_name))) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  stop("Cannot determine script directory for `", script_name, "`. Run from the script folder or with Rscript.", call. = FALSE)
}
script_dir <- resolve_script_dir("Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_pdf <- file.path(output_dir, "Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation.pdf")

stopifnot(file.exists(input_stable6))
stopifnot(file.exists(input_clinical))
stopifnot(file.exists(input_proteome))
stopifnot(file.exists(input_glyco))
if (!clinical_sheet %in% readxl::excel_sheets(input_clinical)) {
  stop("Clinical sheet not found in STable1.xlsx: ", clinical_sheet, call. = FALSE)
}

target_glycopeptide <- "ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S1G0"
target_protein <- "CD274"

zscore <- function(x) {
  as.numeric(scale(as.numeric(x)))
}

extract_feature_row <- function(data, feature_column, feature_id, sample_ids) {
  feature_index <- match(feature_id, data[[feature_column]])
  
  if (is.na(feature_index)) {
    stop(
      "Feature not found: ",
      feature_id,
      " in column ",
      feature_column,
      call. = FALSE
    )
  }
  
  available_ids <- intersect(sample_ids, colnames(data))
  
  if (length(available_ids) == 0) {
    stop(
      "No requested sample IDs were found in the data matrix for feature: ",
      feature_id,
      call. = FALSE
    )
  }
  
  values <- as.numeric(data[feature_index, available_ids, drop = TRUE])
  names(values) <- available_ids
  
  values
}

subtypes <- readxl::read_excel(input_stable6, sheet = "Subtype-cDisc")
subtypes <- as.data.frame(subtypes)
subtypes <- subtypes[, c("id", "protein.subtype")]
colnames(subtypes) <- c("id", "subtype")

subtypes$subtype <- factor(
  subtypes$subtype,
  levels = c("M1", "F1", "M2", "F2", "M3", "F3")
)

clinical <- readxl::read_excel(
  input_clinical,
  sheet = clinical_sheet
)
clinical <- as.data.frame(clinical, check.names = FALSE)

clinical <- clinical[, c("id", "cDisc_age", "cDisc_Gender")]
colnames(clinical) <- c("id", "age", "sex_clinical")

plot_data <- merge(subtypes, clinical, by = "id", all.x = TRUE, sort = FALSE)

plot_data$sex <- ifelse(grepl("^M", plot_data$subtype), "Male", "Female")
plot_data$sex <- ifelse(is.na(plot_data$sex), plot_data$sex_clinical, plot_data$sex)
plot_data$sex <- factor(plot_data$sex, levels = c("Female", "Male"))

proteome <- read.delim(
  input_proteome,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

glyco <- read.delim(
  input_glyco,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_ids <- plot_data$id

cd274_values <- extract_feature_row(
  proteome,
  "ApprovedGeneSymbol",
  target_protein,
  sample_ids
)

itgav_values <- extract_feature_row(
  glyco,
  "Gene.Sequence",
  target_glycopeptide,
  sample_ids
)

common_ids <- intersect(names(cd274_values), names(itgav_values))

measurement_data <- data.frame(
  id = common_ids,
  cd274_raw = cd274_values[common_ids],
  itgav_sial_raw = itgav_values[common_ids],
  stringsAsFactors = FALSE
)

measurement_data$cd274 <- zscore(measurement_data$cd274_raw)
measurement_data$itgav.sial.v2 <- zscore(measurement_data$itgav_sial_raw)

plot_data <- merge(plot_data, measurement_data, by = "id", all.x = FALSE, sort = FALSE)

plot_data <- plot_data[!is.na(plot_data$age) & plot_data$age < 40, ]
plot_data <- plot_data[
  !is.na(plot_data$cd274) &
    !is.na(plot_data$itgav.sial.v2) &
    !is.na(plot_data$sex),
]

if (nrow(plot_data) < 3) {
  stop("Too few complete samples to draw Figure 6N.", call. = FALSE)
}

sex_col <- c(
  "Female" = "#CC0000",
  "Male" = "#0000CC"
)

p <- ggplot(
  plot_data,
  aes(
    x = itgav.sial.v2,
    y = cd274,
    color = sex,
    shape = sex,
    lty = sex
  )
) +
  geom_point() +
  geom_smooth(method = "lm", linetype = 1, fill = NA) +
  stat_cor(
    aes(group = 1),
    label.x.npc = "left",
    label.y.npc = "bottom",
    p.accuracy = 0.001,
    r.accuracy = 0.01
  ) +
  stat_cor(
    label.x.npc = "left",
    label.y.npc = "top",
    p.accuracy = 0.001,
    r.accuracy = 0.01
  ) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    color = "black",
    linetype = 2,
    size = 1,
    se = FALSE
  ) +
  scale_color_manual(values = sex_col) +
  scale_shape_manual(values = c("Female" = 8, "Male" = 18)) +
  theme_bw() +
  ylab("CD274/PD-L1 Protein") +
  xlab("ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S1G0 (Glyco)") +
  theme(legend.position = "right")

pdf(output_pdf, height = 6, width = 6, useDingbats = FALSE)
print(p)
dev.off()

message("Wrote Figure 6N panel to: ", output_pdf)
message(
  "Samples plotted: ",
  nrow(plot_data),
  " (Female = ",
  sum(plot_data$sex == "Female"),
  ", Male = ",
  sum(plot_data$sex == "Male"),
  ")"
)
