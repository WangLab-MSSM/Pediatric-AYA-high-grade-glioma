#!/usr/bin/env Rscript

# ==============================================================================
# File: figure6m_s6o_s6j_6i_itgav_pla2g4a_subtype_analysis.R
# ==============================================================================
#
# Title:
#   Figure 6M, Figure S6O, Figure S6J, and Figure 6I —
#   ITGAV and PLA2G4A multi-omic abundance across sex-stratified
#   proteomic subtypes
#
# Author:
#   Nicole L. Tignor, PhD
#   Department of Genetics and Genomics
#   Icahn School of Medicine at Mount Sinai
#
# Description:
#   This script generates subtype-based visualizations for selected ITGAV and
#   PLA2G4A molecular features across sex-stratified proteomic subtypes in the
#   Discovery cohort. All plotted analyses are restricted to samples with
#   cDisc_age <= 62.
#
#   Figure 6M:
#     - ITGAV sialylated glycopeptide abundance
#       ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S1G0
#
#   Figure S6O:
#     - ITGAV RNA abundance
#     - ITGAV protein abundance
#     - ITGAV unsialylated glycopeptide abundance
#       ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S0G0
#     - ITGAV sialylated glycopeptide abundance
#       ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S1G0
#
#   Figure S6J:
#     - PLA2G4A RNA abundance
#
#   Figure 6I:
#     - PLA2G4A protein abundance
#     - PLA2G4A S377 phosphosite abundance
#
# Clinical metadata mapping:
#   original hclini.m / HARMONY_age
#     -> data/STable1.xlsx / ClinicalTable / cDisc_age
#
# Input files:
#   - data/STable1.xlsx
#   - data/STable6.xlsx
#   - data/cDisc_rna_coding_10192023.tsv
#   - data/cDisc_proteome_imputed_data_09152023.tsv
#   - data/cDisc_phosphosite_imputed_data_ischemia_removed_motif_11032023.tsv
#   - data/Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv
#
# Output files:
#   - figures/Figure6M_ITGAV_sialylated_glycopeptide.pdf
#   - figures/FigureS6O_ITGAV_multiomic_subtype_replication.pdf
#   - figures/FigureS6J_PLA2G4A_RNA_subtype_replication.pdf
#   - figures/Figure6I_PLA2G4A_protein_phosphosite_subtype.pdf
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(readr)
  library(ggplot2)
  library(ggpubr)
})

# ------------------------------------------------------------------------------
# Files and parameters
# ------------------------------------------------------------------------------

clinical_file <- "data/STable1.xlsx"
clinical_sheet <- "ClinicalTable"

subtype_file <- "data/STable6.xlsx"
subtype_sheet <- "Subtype-cDisc"

rna_file <- "data/cDisc_rna_coding_10192023.tsv"
protein_file <- "data/cDisc_proteome_imputed_data_09152023.tsv"
phosphosite_file <- "data/cDisc_phosphosite_imputed_data_ischemia_removed_motif_11032023.tsv"
glyco_file <- "data/Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv"

output_figure6m_pdf <- "figures/Figure6M_ITGAV_sialylated_glycopeptide.pdf"
output_figure_s6o_pdf <- "figures/FigureS6O_ITGAV_multiomic_subtype_replication.pdf"
output_figure_s6j_pdf <- "figures/FigureS6J_PLA2G4A_RNA_subtype_replication.pdf"
output_figure6i_pdf <- "figures/Figure6I_PLA2G4A_protein_phosphosite_subtype.pdf"

age_upper_limit <- 62

itgav_gene <- "ITGAV"
pla2g4a_gene <- "PLA2G4A"
pla2g4a_phosphosite <- "S377"

itgav_glycopeptides <- c(
  "ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S0G0",
  "ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S1G0"
)

subtype_levels <- c("M1", "F1", "M2", "F2", "M3", "F3")

cluster_col <- c(
  "1" = "#B94D3C",
  "2" = "#65B023",
  "3" = "#855F82"
)

subtype_col <- c(
  "M1" = cluster_col[["1"]],
  "F1" = cluster_col[["1"]],
  "M2" = cluster_col[["2"]],
  "F2" = cluster_col[["2"]],
  "M3" = cluster_col[["3"]],
  "F3" = cluster_col[["3"]]
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

sample_id_from_col <- function(x) {
  gsub("\\.", "-", gsub("^((X|A.|G.|P.))", "", x))
}

z_score <- function(x) {
  as.numeric(scale(as.numeric(x)))
}

get_feature_values <- function(data, feature_id, feature_col, sample_ids) {
  row_idx <- which(data[[feature_col]] == feature_id)
  
  if (length(row_idx) != 1) {
    stop(
      "Expected exactly one row for ",
      feature_id,
      " in ",
      feature_col,
      "; found ",
      length(row_idx),
      call. = FALSE
    )
  }
  
  annotation_cols <- c(
    "Modification",
    "ApprovedGeneSymbol",
    "Symbol.V5",
    "OldSymbol",
    "Gene",
    "Sequence",
    "Gene.Sequence",
    "GlycanType",
    "gene.type",
    "coding_gene",
    "Observation_Rate_HOPE",
    "Site",
    "Peptide_GBM",
    "Peptide_HOPE",
    "Motif_GBM",
    "Motif_HOPE"
  )
  
  sample_cols <- setdiff(colnames(data), annotation_cols)
  converted_ids <- sample_id_from_col(sample_cols)
  matched_cols <- sample_cols[match(sample_ids, converted_ids)]
  
  values <- rep(NA_real_, length(sample_ids))
  names(values) <- sample_ids
  
  matched <- !is.na(matched_cols)
  values[matched] <- as.numeric(data[row_idx, matched_cols[matched]])
  
  values
}

get_phosphosite_values <- function(data, gene, site, sample_ids) {
  row_idx <- which(
    data[["ApprovedGeneSymbol"]] == gene &
      grepl(paste0("_", site, "$"), data[["Site"]])
  )
  
  if (length(row_idx) != 1) {
    stop(
      "Expected exactly one row for ",
      gene,
      "_",
      site,
      " in phosphosite data; found ",
      length(row_idx),
      call. = FALSE
    )
  }
  
  get_feature_values(
    data = data,
    feature_id = data[["Site"]][row_idx],
    feature_col = "Site",
    sample_ids = sample_ids
  )
}

make_subtype_boxplot <- function(plot_data) {
  ggplot(
    plot_data,
    aes(x = subtype, y = value, fill = subtype)
  ) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = subtype_col, drop = FALSE) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.8, seed = 123),
      size = 0.5
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_bw() +
    facet_wrap(~ variable, nrow = 1) +
    stat_compare_means(
      ref.group = ".all.",
      label = "p.signif",
      hide.ns = TRUE,
      label.y.npc = "bottom",
      method = "t.test",
      paired = FALSE,
      hjust = 0.5,
      color = "black",
      size = 5,
      fontface = "bold"
    ) +
    stat_compare_means(
      comparisons = list(c("M2", "F2")),
      label = "p.signif",
      hide.ns = TRUE,
      method = "t.test",
      paired = FALSE,
      color = "black",
      size = 5,
      fontface = "bold"
    ) +
    stat_compare_means(
      comparisons = list(c("M1", "F1")),
      label = "p.signif",
      hide.ns = TRUE,
      method = "t.test",
      paired = FALSE,
      color = "black",
      size = 5,
      fontface = "bold"
    ) +
    stat_compare_means(
      comparisons = list(c("M3", "F3")),
      label = "p.signif",
      hide.ns = TRUE,
      method = "t.test",
      paired = FALSE,
      color = "black",
      size = 5,
      fontface = "bold"
    ) +
    labs(
      x = NULL,
      y = "Z-score",
      fill = NULL
    ) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 8),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.margin = margin(t = 15, r = 5, b = 5, l = 5)
    )
}

save_pdf <- function(plot, output_pdf, width, height) {
  dir.create(dirname(output_pdf), showWarnings = FALSE, recursive = TRUE)
  
  pdf(output_pdf, width = width, height = height, useDingbats = FALSE)
  print(plot)
  dev.off()
  
  message("Wrote: ", output_pdf)
}

# ------------------------------------------------------------------------------
# Read metadata
# ------------------------------------------------------------------------------

clinical <- as.data.frame(read_excel(clinical_file, sheet = clinical_sheet))
subtypes <- as.data.frame(read_excel(subtype_file, sheet = subtype_sheet))

clinical$id <- as.character(clinical$id)
clinical$cDisc_age <- as.numeric(clinical$cDisc_age)

subtypes$id <- as.character(subtypes$id)
subtypes$subtype <- factor(subtypes$protein.subtype, levels = subtype_levels)

metadata <- clinical %>%
  transmute(
    id = as.character(id),
    age = cDisc_age
  ) %>%
  inner_join(
    subtypes %>%
      transmute(
        id = as.character(id),
        subtype = factor(protein.subtype, levels = subtype_levels)
      ),
    by = "id"
  ) %>%
  filter(
    !is.na(age),
    age <= age_upper_limit,
    !is.na(subtype)
  )

sample_ids <- metadata$id

# ------------------------------------------------------------------------------
# Read omics matrices
# ------------------------------------------------------------------------------

rna_data <- read_tsv(rna_file, show_col_types = FALSE)
protein_data <- read_tsv(protein_file, show_col_types = FALSE)
phosphosite_data <- read_tsv(phosphosite_file, show_col_types = FALSE)
glyco_data <- read_tsv(glyco_file, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# Extract and scale ITGAV features
# ------------------------------------------------------------------------------

itgav_plot_data_wide <- metadata %>%
  mutate(
    protein_ITGAV = z_score(
      get_feature_values(
        data = protein_data,
        feature_id = itgav_gene,
        feature_col = "ApprovedGeneSymbol",
        sample_ids = sample_ids
      )
    ),
    
    rna_ITGAV = z_score(
      get_feature_values(
        data = rna_data,
        feature_id = itgav_gene,
        feature_col = "ApprovedGeneSymbol",
        sample_ids = sample_ids
      )
    ),
    
    glyco_ITGAV_ANTTQPGIVEGGQVLK_N3H5F1S0G0 = z_score(
      get_feature_values(
        data = glyco_data,
        feature_id = itgav_glycopeptides[1],
        feature_col = "Gene.Sequence",
        sample_ids = sample_ids
      )
    ),
    
    glyco_ITGAV_ANTTQPGIVEGGQVLK_N3H5F1S1G0 = z_score(
      get_feature_values(
        data = glyco_data,
        feature_id = itgav_glycopeptides[2],
        feature_col = "Gene.Sequence",
        sample_ids = sample_ids
      )
    )
  )

# ------------------------------------------------------------------------------
# Convert ITGAV data to long plotting format
# ------------------------------------------------------------------------------

itgav_plot_data <- itgav_plot_data_wide %>%
  pivot_longer(
    cols = c(
      rna_ITGAV,
      protein_ITGAV,
      glyco_ITGAV_ANTTQPGIVEGGQVLK_N3H5F1S0G0,
      glyco_ITGAV_ANTTQPGIVEGGQVLK_N3H5F1S1G0
    ),
    names_to = "variable",
    values_to = "value"
  ) %>%
  filter(
    !is.na(value),
    !is.na(subtype)
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = c(
        "rna_ITGAV",
        "protein_ITGAV",
        "glyco_ITGAV_ANTTQPGIVEGGQVLK_N3H5F1S0G0",
        "glyco_ITGAV_ANTTQPGIVEGGQVLK_N3H5F1S1G0"
      ),
      labels = c(
        "rna_ITGAV",
        "protein_ITGAV",
        "glyco_ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S0G0",
        "glyco_ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S1G0"
      )
    )
  )

# ------------------------------------------------------------------------------
# Extract and scale PLA2G4A features
# ------------------------------------------------------------------------------

pla2g4a_plot_data_wide <- metadata %>%
  mutate(
    protein_PLA2G4A = z_score(
      get_feature_values(
        data = protein_data,
        feature_id = pla2g4a_gene,
        feature_col = "ApprovedGeneSymbol",
        sample_ids = sample_ids
      )
    ),
    
    phospho_PLA2G4A_S377 = z_score(
      get_phosphosite_values(
        data = phosphosite_data,
        gene = pla2g4a_gene,
        site = pla2g4a_phosphosite,
        sample_ids = sample_ids
      )
    ),
    
    rna_PLA2G4A = z_score(
      get_feature_values(
        data = rna_data,
        feature_id = pla2g4a_gene,
        feature_col = "ApprovedGeneSymbol",
        sample_ids = sample_ids
      )
    )
  )

# ------------------------------------------------------------------------------
# Convert PLA2G4A data to long plotting format
# ------------------------------------------------------------------------------

pla2g4a_plot_data <- pla2g4a_plot_data_wide %>%
  pivot_longer(
    cols = c(
      rna_PLA2G4A,
      protein_PLA2G4A,
      phospho_PLA2G4A_S377
    ),
    names_to = "variable",
    values_to = "value"
  ) %>%
  filter(
    !is.na(value),
    !is.na(subtype)
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = c(
        "rna_PLA2G4A",
        "protein_PLA2G4A",
        "phospho_PLA2G4A_S377"
      ),
      labels = c(
        "rna_PLA2G4A",
        "protein_PLA2G4A",
        "phospho_PLA2G4A_S377"
      )
    )
  )

# ------------------------------------------------------------------------------
# Generate Figure 6M
# ------------------------------------------------------------------------------

figure6m_data <- itgav_plot_data %>%
  filter(variable == "glyco_ITGAV_ANTTQPGIVEGGQVLK-N3H5F1S1G0")

figure6m <- make_subtype_boxplot(figure6m_data)

# ------------------------------------------------------------------------------
# Generate Figure S6O
# ------------------------------------------------------------------------------

figure_s6o <- make_subtype_boxplot(itgav_plot_data)

# ------------------------------------------------------------------------------
# Generate Figure S6J
# ------------------------------------------------------------------------------

figure_s6j_data <- pla2g4a_plot_data %>%
  filter(variable == "rna_PLA2G4A")

figure_s6j <- make_subtype_boxplot(figure_s6j_data)

# ------------------------------------------------------------------------------
# Generate Figure 6I
# ------------------------------------------------------------------------------

figure6i_data <- pla2g4a_plot_data %>%
  filter(variable %in% c("protein_PLA2G4A", "phospho_PLA2G4A_S377")) %>%
  mutate(
    variable = factor(
      variable,
      levels = c("protein_PLA2G4A", "phospho_PLA2G4A_S377")
    )
  )

figure6i <- make_subtype_boxplot(figure6i_data)

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

save_pdf(figure6m, output_figure6m_pdf, width = 3.25, height = 4)
save_pdf(figure_s6o, output_figure_s6o_pdf, width = 8.5, height = 4)
save_pdf(figure_s6j, output_figure_s6j_pdf, width = 3.25, height = 4)
save_pdf(figure6i, output_figure6i_pdf, width = 5.5, height = 4)

message("Age upper limit: ", age_upper_limit)
message("Figure 6M plotted samples: ", length(unique(figure6m_data$id)))
message("Figure S6O plotted samples: ", length(unique(itgav_plot_data$id)))
message("Figure S6J plotted samples: ", length(unique(figure_s6j_data$id)))
message("Figure 6I plotted samples: ", length(unique(figure6i_data$id)))