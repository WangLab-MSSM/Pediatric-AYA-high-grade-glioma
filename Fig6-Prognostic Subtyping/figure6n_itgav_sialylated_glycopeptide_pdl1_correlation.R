#!/usr/bin/env Rscript

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
#   figure6n_itgav_sialylated_glycopeptide_pdl1_correlation.pdf
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

# Resolve paths relative to repository root.
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)

script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  NA_character_
}

repo_root <- if (!is.na(script_path)) dirname(dirname(script_path)) else getwd()

if (
  !file.exists(file.path(repo_root, "data", "STable6.xlsx")) &&
  basename(getwd()) == "code_dump"
) {
  repo_root <- dirname(getwd())
}

input_stable6 <- file.path(repo_root, "data", "STable6.xlsx")
input_clinical <- file.path(repo_root, "data", "clinical_data_04032026.tsv")
input_proteome <- file.path(repo_root, "data", "cDisc_proteome_imputed_data_09152023.tsv")
input_glyco <- file.path(repo_root, "data", "Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv")

output_pdf <- file.path(
  repo_root,
  "figure6n_itgav_sialylated_glycopeptide_pdl1_correlation.pdf"
)

stopifnot(file.exists(input_stable6))
stopifnot(file.exists(input_clinical))
stopifnot(file.exists(input_proteome))
stopifnot(file.exists(input_glyco))

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

clinical <- read.delim(
  input_clinical,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

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