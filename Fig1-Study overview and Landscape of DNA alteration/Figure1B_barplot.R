#!/usr/bin/env Rscript

# ============================================================
# Figure 1B | Subtype Distribution Across Development
# File: figure1b_subtype_barplot.R
#
# Description:
#   Generates a stacked bar plot showing subtype proportions
#   in pediatric (PED) and adolescent/young adult (AYA) tumors
#   from the study data object.
#
# Input:
#   - data/pediatric_aya_hgg_study_data.rds
#
# Output:
#   - Figure1B_barplot_subtype.pdf
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
})

# ---- I/O ----
input_file <- "data/pediatric_aya_hgg_study_data.rds"
output_file <- "Figure1B_barplot_subtype.pdf"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# ---- Load data ----
study.data <- readRDS(input_file)
clinical.table <- study.data$clinical

# ---- Preprocess ----

# Collapse ADO and YA into AYA
clinical.table$age_class <- ifelse(
  clinical.table$cDisc_age_class_name_derived %in% c("ADO", "YA"),
  "AYA",
  as.character(clinical.table$cDisc_age_class_name_derived)
)

# Assign flagged samples to "Flag"
clinical.table$cancer_group <- ifelse(
  clinical.table$flag == 1,
  "Flag",
  as.character(clinical.table$Disc_Cancer_Group)
)

# Restrict to PED and AYA
plot.df <- subset(
  clinical.table,
  age_class %in% c("PED", "AYA") & !is.na(cancer_group)
)

# ---- Compute proportions ----
count_mat <- table(plot.df$age_class, plot.df$cancer_group)
prop_mat <- prop.table(t(count_mat), margin = 2)

# Ensure both PED and AYA exist (robustness)
for (grp in c("PED", "AYA")) {
  if (!grp %in% colnames(prop_mat)) {
    prop_mat <- cbind(prop_mat, setNames(data.frame(rep(0, nrow(prop_mat))), grp))
  }
}
prop_mat <- prop_mat[, c("PED", "AYA"), drop = FALSE]

# ---- Format for plotting ----
plot.data <- data.frame(
  subtype = rownames(prop_mat),
  PED = prop_mat[, "PED"],
  AYA = prop_mat[, "AYA"],
  row.names = NULL,
  stringsAsFactors = FALSE
)

plot.data <- melt(
  plot.data,
  id.vars = "subtype",
  variable.name = "age_class",
  value.name = "prop"
)

plot.data$age_class <- factor(plot.data$age_class, levels = c("PED", "AYA"))

plot.data$label <- ifelse(
  plot.data$prop > 0,
  paste0(round(plot.data$prop * 100), "%"),
  ""
)

# ---- Colors ----
color_palette <- c(
  "(DHG) Diffuse Hemispheric Glioma" = "#A5772B",
  "(DMG) Diffuse Midline Glioma" = "#1B3360",
  "(HGG) High Grade Glioma (not otherwise specified)" = "#4393C3",
  "(IHG) Infantile Hemispheric Glioma" = "#E62A8A",
  "(PXA) Pleomorphic Xanthoastrocytoma" = "#660821",
  "Flag" = "gray80"
)

# Check for missing colors
missing_colors <- setdiff(unique(plot.data$subtype), names(color_palette))
if (length(missing_colors) > 0) {
  stop("Missing colors for subtype(s): ", paste(missing_colors, collapse = ", "))
}

# ---- Plot ----
p <- ggplot(plot.data, aes(x = age_class, y = prop, fill = subtype)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.1) +
  coord_flip() +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = color_palette) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  labs(
    x = "Age Class",
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

# ---- Save ----
ggsave(
  filename = output_file,
  plot = p,
  width = 11,
  height = 3
)