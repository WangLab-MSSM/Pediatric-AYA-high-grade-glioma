#!/usr/bin/env Rscript

# ==============================================================================
# Author: Nicole L. Tignor
# Purpose: Recreate Figure 4B from public supplementary table data
# Input: data/STable4.xlsx, sheet SA-Glyco-Disc
# Output: Figure 4B PDF
# ==============================================================================

required_packages <- c(
  "readxl",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install required package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

input_file <- file.path("../data", "STable4.xlsx")
sheet_name <- "SA-Glyco-Disc"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, "Figure4B_survival_assoc_glyco.pdf")
fdr_threshold <- 0.10
signed_fdr_threshold <- -log10(fdr_threshold)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

available_sheets <- readxl::excel_sheets(input_file)
if (!sheet_name %in% available_sheets) {
  stop(
    "Sheet not found: ", sheet_name,
    "\nAvailable sheets: ", paste(available_sheets, collapse = ", "),
    call. = FALSE
  )
}

raw_data <- readxl::read_excel(input_file, sheet = sheet_name, col_types = "text")

required_columns <- c(
  "glycoppetide",
  "membrane.loc",
  "PED.comb.signed.p.female.adj",
  "ADO.comb.signed.p.female.adj",
  "PED.comb.signed.fdr.female.adj",
  "ADO.comb.signed.fdr.female.adj",
  "PED.comb.signed.pmale.adj",
  "ADO.comb.signed.pmale.adj",
  "PED.comb.signed.fdrmale.adj",
  "ADO.comb.signed.fdrmale.adj",
  "PED.comb.signed.p.female.noadj",
  "ADO.comb.signed.p.female.noadj",
  "PED.comb.signed.fdr.female.noadj",
  "ADO.comb.signed.fdr.female.noadj",
  "PED.comb.signed.pmale.noadj",
  "ADO.comb.signed.pmale.noadj",
  "PED.comb.signed.fdrmale.noadj",
  "ADO.comb.signed.fdrmale.noadj"
)

missing_columns <- setdiff(required_columns, names(raw_data))
if (length(missing_columns) > 0) {
  stop(
    "Required column(s) missing from ", sheet_name, ": ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

numeric_columns <- setdiff(required_columns, c("glycoppetide", "membrane.loc"))
for (column in numeric_columns) {
  raw_data[[column]] <- suppressWarnings(as.numeric(raw_data[[column]]))
}

raw_data$gene <- sub("_.*", "", raw_data$glycoppetide)
raw_data$glycan <- sub(".*-", "", raw_data$glycoppetide)
raw_data$display_label <- paste(raw_data$gene, raw_data$glycan, sep = " - ")

published_rows <- data.frame(
  group = c(
    rep("Regulation of\nNervous System\nDevelopment", 11),
    rep("TMP Group MF-24", 3),
    rep("Extracellular\nMatrix\nOrganization", 16),
    rep("TMP Group MF-13", 3),
    rep("Coagulation", 2)
  ),
  display_label = c(
    "EPHB3 - N2H5F0S0G0",
    "L1CAM - N5H3F1S0G0",
    "MEGF8 - N2H8F0S0G0",
    "NPTN - N5H3F0S0G0",
    "NTRK2 - N2H5F0S0G0",
    "PLXNA4 - N2H6F0S0G0",
    "PLXNB2 - N2H8F0S0G0",
    "PLXNB3 - N2H5F0S0G0",
    "PLXNC1 - N2H5F0S0G0",
    "PTPRZ1 - N4H4F0S1G0",
    "THY1 - N3H5F1S0G0",
    "CNTN1 - N2H8F0S0G0",
    "OMG - N2H8F0S0G0",
    "OMG - N2H7F0S0G0",
    "ADAM17 - N2H5F0S0G0",
    "BCAN - N5H3F3S0G0",
    "CD44 - N4H5F0S1G0",
    "CD47 - N2H5F0S0G0",
    "COL6A2 - N4H5F0S1G0",
    "ITGA1 - N2H5F0S0G0",
    "ITGA3 - N2H5F0S0G0",
    "ITGAV - N2H9F0S0G0",
    "ITGB1 - N2H5F0S0G0",
    "ITGB3 - N2H5F0S0G0",
    "LAMC1 - N3H4F0S0G0",
    "NCAM1 - N5H3F1S0G0",
    "NID2 - N4H4F1S1G0",
    "PLOD1 - N2H7F0S0G0",
    "PLOD3 - N2H5F0S0G0",
    "PXDN - N2H5F0S0G0",
    "CD276 - N2H5F0S0G0",
    "LAMP1 - N3H4F1S1G0",
    "LAMP1 - N3H6F1S0G0",
    "LAMP2 - N2H6F0S0G0",
    "LRP1 - N2H8F0S0G0"
  ),
  stringsAsFactors = FALSE
)

fdr_columns <- c(
  "PED.comb.signed.fdrmale.adj",
  "ADO.comb.signed.fdrmale.adj",
  "PED.comb.signed.fdr.female.adj",
  "ADO.comb.signed.fdr.female.adj",
  "PED.comb.signed.fdrmale.noadj",
  "ADO.comb.signed.fdrmale.noadj",
  "PED.comb.signed.fdr.female.noadj",
  "ADO.comb.signed.fdr.female.noadj"
)

p_columns <- c(
  "PED.comb.signed.pmale.adj",
  "ADO.comb.signed.pmale.adj",
  "PED.comb.signed.p.female.adj",
  "ADO.comb.signed.p.female.adj",
  "PED.comb.signed.pmale.noadj",
  "ADO.comb.signed.pmale.noadj",
  "PED.comb.signed.p.female.noadj",
  "ADO.comb.signed.p.female.noadj"
)

resolve_row <- function(label) {
  candidate_index <- which(raw_data$display_label == label)
  if (length(candidate_index) == 0) {
    stop("Published Figure 4B row not found in sheet: ", label, call. = FALSE)
  }

  if (length(candidate_index) == 1) {
    return(candidate_index)
  }

  fdr_score <- apply(
    abs(as.data.frame(raw_data[candidate_index, fdr_columns, drop = FALSE])),
    1,
    function(x) max(x, na.rm = TRUE)
  )
  fdr_score[is.infinite(fdr_score)] <- NA_real_

  if (all(is.na(fdr_score))) {
    p_score <- apply(
      abs(as.data.frame(raw_data[candidate_index, p_columns, drop = FALSE])),
      1,
      function(x) max(x, na.rm = TRUE)
    )
    p_score[is.infinite(p_score)] <- NA_real_
    chosen <- candidate_index[which.max(p_score)]
  } else {
    chosen <- candidate_index[which.max(fdr_score)]
  }

  message(
    "Resolved duplicate printed label '", label, "' by selecting full glycopeptide: ",
    raw_data$glycoppetide[chosen]
  )

  chosen
}

selected_index <- vapply(published_rows$display_label, resolve_row, integer(1))
plot_rows <- raw_data[selected_index, , drop = FALSE]
plot_rows$group <- factor(published_rows$group, levels = unique(published_rows$group))
plot_rows$published_label <- published_rows$display_label

value_columns <- c(
  "PED.comb.signed.pmale.adj",
  "ADO.comb.signed.pmale.adj",
  "PED.comb.signed.p.female.adj",
  "ADO.comb.signed.p.female.adj",
  "PED.comb.signed.pmale.noadj",
  "ADO.comb.signed.pmale.noadj",
  "PED.comb.signed.p.female.noadj",
  "ADO.comb.signed.p.female.noadj"
)

significance_columns <- c(
  "PED.comb.signed.fdrmale.adj",
  "ADO.comb.signed.fdrmale.adj",
  "PED.comb.signed.fdr.female.adj",
  "ADO.comb.signed.fdr.female.adj",
  "PED.comb.signed.fdrmale.noadj",
  "ADO.comb.signed.fdrmale.noadj",
  "PED.comb.signed.fdr.female.noadj",
  "ADO.comb.signed.fdr.female.noadj"
)

heatmap_matrix <- as.matrix(plot_rows[, value_columns, drop = FALSE])
significance_matrix <- as.matrix(plot_rows[, significance_columns, drop = FALSE])
significance_matrix <- ifelse(
  abs(significance_matrix) >= signed_fdr_threshold,
  significance_matrix,
  NA_real_
)

column_labels <- c("PED", "ADO", "PED", "ADO", "PED", "ADO", "PED", "ADO")
colnames(heatmap_matrix) <- column_labels
colnames(significance_matrix) <- column_labels
rownames(heatmap_matrix) <- plot_rows$published_label
rownames(significance_matrix) <- plot_rows$published_label

column_metadata <- data.frame(
  Adjustment = factor(
    c(rep("Adj.", 4), rep("No Adj.", 4)),
    levels = c("Adj.", "No Adj.")
  ),
  Sex = factor(
    rep(c("M", "F"), each = 2, times = 2),
    levels = c("M", "F")
  ),
  `Age class` = factor(
    rep(c("PED", "ADO"), times = 4),
    levels = c("PED", "ADO")
  ),
  check.names = FALSE
)

column_split <- factor(
  paste(column_metadata$Adjustment, column_metadata$Sex),
  levels = unique(paste(column_metadata$Adjustment, column_metadata$Sex))
)

plot_rows$sialylated <- grepl("S[1-9]", plot_rows$glycan)
plot_rows$fucosylated <- grepl("F[1-9]", plot_rows$glycan)
plot_rows$high_mannose <- grepl("^N2H[5-9]F0S0G0$", plot_rows$glycan)
plot_rows$membrane <- plot_rows$membrane.loc == "Yes"

pathway_colors <- c(
  "Regulation of\nNervous System\nDevelopment" = "#b7d900",
  "TMP Group MF-24" = "#c018d9",
  "Extracellular\nMatrix\nOrganization" = "#2ecc71",
  "TMP Group MF-13" = "#2146d0",
  "Coagulation" = "#d7191c"
)

sex_colors <- c("M" = "#0000cc", "F" = "#d00000")
age_class_colors <- c("PED" = "#bfe7bf", "ADO" = "#6fbd6a")
adjustment_colors <- c("Adj." = "#f7ad0b", "No Adj." = "#8e4d07")
binary_colors <- c("TRUE" = "black", "FALSE" = "#eeeeee")

heatmap_colors <- circlize::colorRamp2(
  c(-5, 0, 5),
  rev(colorRampPalette(RColorBrewer::brewer.pal(5, "PiYG"))(3))
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = column_metadata,
  col = list(
    Adjustment = adjustment_colors,
    Sex = sex_colors,
    `Age class` = age_class_colors
  ),
  annotation_name_side = "left",
  show_legend = TRUE
)

left_annotation <- ComplexHeatmap::rowAnnotation(
  Pathway = plot_rows$group,
  Sial = as.character(plot_rows$sialylated),
  Fuco = as.character(plot_rows$fucosylated),
  HM = as.character(plot_rows$high_mannose),
  Mem. = as.character(plot_rows$membrane),
  col = list(
    Pathway = pathway_colors,
    Sial = binary_colors,
    Fuco = binary_colors,
    HM = binary_colors,
    Mem. = c("TRUE" = "#18b957", "FALSE" = "#eeeeee")
  ),
  annotation_name_rot = 90,
  annotation_name_gp = grid::gpar(fontsize = 7),
  simple_anno_size = grid::unit(3.1, "mm"),
  show_legend = TRUE
)

glyco_heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Glyco-peptide\nSurv. Assoc.",
  col = heatmap_colors,
  na_col = "grey88",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = plot_rows$group,
  column_split = column_split,
  row_gap = grid::unit(c(1.5, 1.5, 1.5, 1.5), "mm"),
  column_gap = grid::unit(c(1.2, 2.6, 1.2), "mm"),
  border = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 0.6),
  top_annotation = top_annotation,
  left_annotation = left_annotation,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 6.5),
  row_title_side = "left",
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 7),
  column_names_rot = 90,
  column_title_rot = 45,
  width = grid::unit(ncol(heatmap_matrix) * 4.8, "mm"),
  height = grid::unit(nrow(heatmap_matrix) * 4.1, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(significance_matrix[i, j])) {
      grid::grid.points(
        x,
        y,
        pch = 16,
        size = grid::unit(1.05, "mm"),
        gp = grid::gpar(col = "black")
      )
    }
  }
)

fdr_point_legend <- ComplexHeatmap::Legend(
  title = "Significance",
  labels = paste0("FDR < ", fdr_threshold),
  type = "points",
  pch = 16,
  size = grid::unit(1.05, "mm"),
  legend_gp = grid::gpar(col = "black"),
  title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 7)
)

message(
  "Using signed p-value columns for tile color and signed FDR columns for ",
  "FDR < ", fdr_threshold, " point overlays."
)
message(
  "Figure row grouping/order follows the published panel because pathway/module ",
  "membership is not encoded in SA-Glyco-Disc; all plotted association values ",
  "are read from ", input_file, "."
)

grDevices::cairo_pdf(output_file, width = 7.0, height = 8.2)
grid::grid.newpage()
grid::grid.text(
  "Glyco-peptide Survival Assoc.",
  x = grid::unit(0.57, "npc"),
  y = grid::unit(0.975, "npc"),
  gp = grid::gpar(fontsize = 13, fontface = "bold")
)
ComplexHeatmap::draw(
  glyco_heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(fdr_point_legend),
  newpage = FALSE,
  padding = grid::unit(c(10, 4, 4, 4), "mm")
)
grDevices::dev.off()

message("Saved ", output_file)
