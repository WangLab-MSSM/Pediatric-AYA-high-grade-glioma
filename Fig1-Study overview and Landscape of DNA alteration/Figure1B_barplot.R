#!/usr/bin/env Rscript

# ==============================================================================
# File: figure1b_subtype_barplot.R
# ==============================================================================
#
# Title:
#   Figure 1B — Subtype distribution across developmental age groups
#
# Author:
#   Nicole L. Tignor, PhD
#   Department of Genetics and Genomics
#   Icahn School of Medicine at Mount Sinai
#
# Description:
#   Generates a stacked bar plot showing tumor subtype proportions in pediatric
#   (PED) and adolescent/young adult (AYA) tumors using the clinical data file
#   directly from the data directory. 
#
# Input file:
#   - data/cDisc_clinical_data_04032026.tsv
#
# Output file:
#   - Figure1B_barplot_subtype.pdf
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------------------------
# Input / output
# ------------------------------------------------------------------------------

clinical_file <- "../data/cDisc_clinical_data_04032026.tsv"
output_file <- "Figure1B_barplot_subtype.pdf"

if (!file.exists(clinical_file)) {
  stop("Clinical file not found: ", clinical_file)
}

# ------------------------------------------------------------------------------
# Read clinical data
# ------------------------------------------------------------------------------

clinical_table <- read_tsv(
  clinical_file,
  show_col_types = FALSE
) %>%
  as.data.frame()

# ------------------------------------------------------------------------------
# Preprocess
# ------------------------------------------------------------------------------

clinical_table <- clinical_table %>%
  mutate(
    age_class = ifelse(
      cDisc_age_class_name_derived %in% c("ADO", "YA"),
      "AYA",
      as.character(cDisc_age_class_name_derived)
    ),
    cancer_group = ifelse(
      flag == 1,
      "Flag",
      as.character(Disc_Cancer_Group)
    ),
    cancer_group = ifelse(
      cancer_group %in% c("NA", "", "NaN"),
      NA,
      cancer_group
    )
  )

plot_df <- clinical_table %>%
  filter(
    age_class %in% c("PED", "AYA"),
    !is.na(cancer_group)
  )

# ------------------------------------------------------------------------------
# Compute proportions
# ------------------------------------------------------------------------------

plot_data <- plot_df %>%
  count(age_class, cancer_group, name = "n") %>%
  group_by(age_class) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup() %>%
  complete(
    age_class = c("PED", "AYA"),
    cancer_group,
    fill = list(n = 0, prop = 0)
  ) %>%
  mutate(
    age_class = factor(age_class, levels = c("PED", "AYA")),
    label = ifelse(prop > 0, paste0(round(prop * 100), "%"), "")
  )

# ------------------------------------------------------------------------------
# Colors
# ------------------------------------------------------------------------------

color_palette <- c(
  "(DHG) Diffuse Hemispheric Glioma" = "#A5772B",
  "(DMG) Diffuse Midline Glioma" = "#1B3360",
  "(HGG) High Grade Glioma (not otherwise specified)" = "#4393C3",
  "(IHG) Infantile Hemispheric Glioma" = "#E62A8A",
  "(PXA) Pleomorphic Xanthoastrocytoma" = "#660821",
  "Flag" = "gray80"
)

missing_colors <- setdiff(unique(plot_data$cancer_group), names(color_palette))

if (length(missing_colors) > 0) {
  stop(
    "Missing colors for subtype(s): ",
    paste(missing_colors, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = age_class,
    y = prop,
    fill = cancer_group
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    color = "black",
    linewidth = 0.1
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = percent_format(),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = color_palette) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  labs(
    x = "Age class",
    y = "Proportion",
    fill = "Subtype"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

ggsave(
  filename = output_file,
  plot = p,
  width = 11,
  height = 3,
  useDingbats = FALSE
)

message("Wrote: ", output_file)