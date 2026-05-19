#!/usr/bin/env Rscript
# Author: Nicole L. Tignor
# Purpose: Recreate Figure 2A tumor/normal AD-TMP heatmap from public data using ageTMP.
#
# Design note:
# ------------
# This script generates the Figure 2A tumor-normal AD-TMP heatmap from
# public study data and package-provided reference data. The analysis places
# tumor age-dependent molecular trajectories and normal developmental
# trajectories in a shared protein feature space while preserving the
# trajectory structure inferred separately within each biological context.
#
# Tumor and normal AD-TMPs are estimated independently, clustered within their
# respective contexts, and aligned by matched protein/sex rows. This design
# allows the heatmap to show how age-associated tumor programs correspond to,
# diverge from, or reorganize normal developmental molecular programs.
#
# Curated manuscript annotations, including highlighted proteins and pathway
# membership flags, are applied at the visualization stage. The analytic
# matrices, trajectory estimates, and cluster assignments are generated from
# documented inputs through the ageTMP workflow.

required <- c("ComplexHeatmap", "circlize", "grid", "readxl", "RColorBrewer")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("ageTMP")) {
  pkgload::load_all("ageTMP", quiet = TRUE)
} else if (requireNamespace("ageTMP", quietly = TRUE)) {
  library(ageTMP)
} else {
  stop("Install ageTMP or run this script from the repository root containing ageTMP/.", call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

repo_data <- "../data"
clinical_file <- file.path(repo_data, "clinical_data_04032026.tsv")
protein_file <- file.path(repo_data, "cDisc_proteome_imputed_data_09152023.tsv")
stable2_file <- file.path(repo_data, "STable2.xlsx")
for (file in c(clinical_file, protein_file, stable2_file)) {
  if (!file.exists(file)) {
    stop("Required input file not found: ", file, call. = FALSE)
  }
}

stable2_sheets <- readxl::excel_sheets(stable2_file)
for (sheet in c("TTMP_Protein_Group", "TN_Diff")) {
  if (!sheet %in% stable2_sheets) {
    stop("Required STable2 sheet not found: ", sheet, call. = FALSE)
  }
}

message("Loading public clinical, protein, and STable2 data...")
clinical <- ageTMP_load_clinical(repo_data)
protein_raw <- ageTMP_load_molecular(repo_data, "protein")
groups <- as.data.frame(ageTMP_load_supplement(repo_data, "STable2", "TTMP_Protein_Group"))
tn_diff <- as.data.frame(ageTMP_load_supplement(repo_data, "STable2", "TN_Diff"))

required_clinical <- c("id", "cDisc_age", "cDisc_Gender", "cDisc_age_class_name_derived")
missing_clinical <- setdiff(required_clinical, names(clinical))
if (length(missing_clinical) > 0) {
  stop("Clinical data missing required column(s): ", paste(missing_clinical, collapse = ", "), call. = FALSE)
}
required_group <- c("Gene", "Male T-TMP Group", "Female T-TMP Group", "Male:Female T-TMP Group")
missing_group <- setdiff(required_group, names(groups))
if (length(missing_group) > 0) {
  stop("TTMP_Protein_Group missing required column(s): ", paste(missing_group, collapse = ", "), call. = FALSE)
}
required_protein <- c("ApprovedGeneSymbol", "Symbol.V5")
missing_protein <- setdiff(required_protein, names(protein_raw))
if (length(missing_protein) > 0) {
  stop("Protein data missing required annotation column(s): ", paste(missing_protein, collapse = ", "), call. = FALSE)
}

message("Preparing gene-level protein matrix...")
annotation_cols <- match(c("ApprovedGeneSymbol", "Symbol.V5", "OldSymbol", "coding_gene"), names(protein_raw))
annotation_cols <- annotation_cols[!is.na(annotation_cols)]
protein_parts <- ageTMP_split_annotation_matrix(protein_raw, annotation_cols = annotation_cols)
protein_mat <- ageTMP_collapse_matrix_by_feature(
  protein_parts$matrix,
  protein_parts$annotation$ApprovedGeneSymbol
)

clinical$id <- ageTMP_normalize_sample_ids(clinical$id)
clinical$age <- as.numeric(clinical$cDisc_age)
clinical$sex <- clinical$cDisc_Gender
clinical$age_class <- clinical$cDisc_age_class_name_derived

sample_ids <- intersect(colnames(protein_mat), clinical$id)
protein_mat <- protein_mat[, sample_ids, drop = FALSE]
clinical <- clinical[match(sample_ids, clinical$id), , drop = FALSE]

age_levels <- c("PED", "ADO", "YA", "ADULT")
plot_clinical <- clinical[
  !is.na(clinical$age) &
    !is.na(clinical$age_class) &
    clinical$age_class %in% age_levels &
    clinical$id %in% colnames(protein_mat),
  , drop = FALSE
]
plot_clinical$age_class <- factor(plot_clinical$age_class, levels = age_levels)
plot_clinical <- plot_clinical[order(plot_clinical$age), , drop = FALSE]

features <- intersect(groups$Gene, rownames(protein_mat))
groups <- groups[match(features, groups$Gene), , drop = FALSE]
if (length(features) == 0) {
  stop("No STable2 T-TMP genes overlap the public protein matrix.", call. = FALSE)
}
message("Restricting Figure 2A heatmap to proteins present in tumor and normal/reference matrices...")
normal_reference <- ageTMP_load_normal_reference()
normal_protein <- normal_reference$protein$matrix
normal_metadata <- normal_reference$protein$sample_metadata
groups_full <- groups
features_all_tumor <- intersect(groups_full$Gene, rownames(protein_mat))
features <- intersect(features_all_tumor, rownames(normal_protein))
groups <- groups_full[match(features, groups_full$Gene), , drop = FALSE]
if (length(features) == 0) {
  stop("No STable2 T-TMP genes overlap both the public protein matrix and package normal reference.", call. = FALSE)
}
message("Fitting tumor/normal AD-TMP trajectories for ", length(features), " proteins.")
message("Using fixed span = 1.5 for this public-data reconstruction. Published protein_tadj50 used adaptive feature/sex-specific spans.")

prediction_ages <- sort(unique(plot_clinical$age))
prediction_ids <- sprintf("age_%03d", seq_along(prediction_ages))
age_class_from_age <- function(age) {
  out <- rep(NA_character_, length(age))
  out[age < 15] <- "PED"
  out[age >= 15 & age < 26] <- "ADO"
  out[age >= 26 & age < 40] <- "YA"
  out[age >= 40] <- "ADULT"
  factor(out, levels = age_levels)
}
col_age_class <- age_class_from_age(prediction_ages)

trajectory <- ageTMP_compare_normal_tumor_trajectory(
  tumor_mat = protein_mat,
  tumor_metadata = clinical,
  normal_mat = normal_protein,
  normal_metadata = normal_metadata,
  features = features,
  tumor_sample_col = "id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  normal_sample_col = "ID",
  normal_age_col = "Age",
  normal_sex_col = "Gender",
  normal_covariates = c("pH", "PMI", "Ethnicity"),
  center_age_range = c(0, 50),
  span = 1.5,
  prediction_ages = prediction_ages,
  ci_level = 0.95
)
trajectory$sample_id <- prediction_ids[match(trajectory$age, prediction_ages)]

make_fit_matrix <- function(traj, sex, tissue, feature_order, sample_order) {
  sex_df <- traj[traj$sex == sex & traj$tissue == tissue, c("feature", "sample_id", "fit"), drop = FALSE]
  out <- matrix(NA_real_, nrow = length(feature_order), ncol = length(sample_order),
                dimnames = list(feature_order, sample_order))
  split_df <- split(sex_df, sex_df$feature)
  for (feature in names(split_df)) {
    idx <- match(split_df[[feature]]$sample_id, sample_order)
    out[feature, idx] <- split_df[[feature]]$fit
  }
  out
}

cluster_normal_rows <- function(mat, k = 4) {
  scaled <- t(scale(t(mat)))
  scaled[!is.finite(scaled)] <- 0
  set.seed(1)
  km <- stats::kmeans(scaled, centers = k, nstart = 50, iter.max = 100)
  centers <- km$centers
  age_score <- rowMeans(centers[, seq_len(max(1, floor(ncol(centers) / 4))), drop = FALSE]) -
    rowMeans(centers[, seq.int(ceiling(ncol(centers) * 3 / 4), ncol(centers)), drop = FALSE])
  remap <- setNames(seq_len(k), names(sort(age_score, decreasing = TRUE)))
  out <- factor(remap[as.character(km$cluster)], levels = as.character(seq_len(k)))
  names(out) <- rownames(mat)
  out
}

message("Assigning male and female N-TMP trajectory groups...")
normal_cluster_feature_order <- groups$Gene
normal_consensus_k4 <- normal_reference$protein$clusters$consensus_k4 %||% NULL
if (!is.null(normal_consensus_k4)) {
  message("Using packaged manuscript N-TMP consensus clusters with published remapping.")
  normal_cluster_by_sex <- list(
    Male = setNames(
      factor(
        as.character(normal_consensus_k4[match(paste0("Male,", normal_cluster_feature_order), names(normal_consensus_k4))]),
        levels = as.character(1:4)
      ),
      normal_cluster_feature_order
    ),
    Female = setNames(
      factor(
        as.character(normal_consensus_k4[match(paste0("Female,", normal_cluster_feature_order), names(normal_consensus_k4))]),
        levels = as.character(1:4)
      ),
      normal_cluster_feature_order
    )
  )
} else {
  message("Packaged N-TMP consensus clusters not found; estimating joint k-means clusters as fallback.")
  # The fallback clusters are estimated in one shared male+female space. This is
  # the computational counterpart of the Sankey logic in the manuscript: a color
  # such as normal group 2 must have the same trajectory meaning in the male and
  # female normal matrices, rather than being independently discovered and
  # independently colored in each sex.
  male_normal_for_cluster <- make_fit_matrix(
    trajectory, "Male", "Normal", normal_cluster_feature_order, prediction_ids
  )
  female_normal_for_cluster <- make_fit_matrix(
    trajectory, "Female", "Normal", normal_cluster_feature_order, prediction_ids
  )
  joint_normal_for_cluster <- rbind(male_normal_for_cluster, female_normal_for_cluster)
  rownames(joint_normal_for_cluster) <- c(
    paste0("Male,", normal_cluster_feature_order),
    paste0("Female,", normal_cluster_feature_order)
  )
  joint_normal_clusters <- cluster_normal_rows(joint_normal_for_cluster, k = 4)
  normal_cluster_by_sex <- list(
    Male = setNames(
      joint_normal_clusters[match(paste0("Male,", normal_cluster_feature_order), names(joint_normal_clusters))],
      normal_cluster_feature_order
    ),
    Female = setNames(
      joint_normal_clusters[match(paste0("Female,", normal_cluster_feature_order), names(joint_normal_clusters))],
      normal_cluster_feature_order
    )
  )
}
normal_alignment_summary <- table(
  male = normal_cluster_by_sex$Male[normal_cluster_feature_order],
  female = normal_cluster_by_sex$Female[normal_cluster_feature_order]
)
print(normal_alignment_summary)

mark_genes <- c(
  "ALDH9A1", "ATL1", "CACNA2D1", "CALB2", "CNTN1", "CYC1",
  "DCLK1", "EPB41L3", "GNAZ", "L1CAM", "MAPT", "MGST3",
  "NCAM2", "NDUFA6", "NDUFB4", "OPCML", "OXCT1", "PDE2A",
  "SLC8A2", "VDAC1"
)

make_sex_panel <- function(sex) {
  cl_col <- if (sex == "Male") "Male T-TMP Group" else "Female T-TMP Group"
  opposite_col <- if (sex == "Male") "Female T-TMP Group" else "Male T-TMP Group"

  sankey_base <- data.frame(
    gene = groups$Gene,
    male.cl = factor(as.character(groups[["Male T-TMP Group"]]), levels = c("1", "2", "3", "4")),
    female.cl = factor(as.character(groups[["Female T-TMP Group"]]), levels = c("1", "2", "3", "4")),
    n.male.cl = factor(as.character(normal_cluster_by_sex$Male[groups$Gene]), levels = c("1", "2", "3", "4")),
    n.female.cl = factor(as.character(normal_cluster_by_sex$Female[groups$Gene]), levels = c("1", "2", "3", "4")),
    stringsAsFactors = FALSE
  )
  sankey_base$join <- paste0(sankey_base$male.cl, sankey_base$n.male.cl, sankey_base$female.cl, sankey_base$n.female.cl)

  # Sankey-aligned row ordering.
  #
  # This mirrors the manuscript code more literally than a per-sex sort. A
  # single gene order is built from the full tumor/normal tuple:
  #   male T-TMP, male N-TMP, female T-TMP, female N-TMP.
  # That global order is then reused for both sex-specific heatmaps. ComplexHeatmap
  # performs the visible row splitting by the sex-specific T-TMP group, while the
  # within-split order still carries the full Sankey relationship.
  sankey_order <- order(sankey_base$join, seq_len(nrow(sankey_base)))
  sankey_base <- sankey_base[sankey_order, , drop = FALSE]
  tumor_features <- sankey_base$gene

  ncl <- if (sex == "Male") sankey_base$n.male.cl else sankey_base$n.female.cl
  tumor_base <- data.frame(
    gene = sankey_base$gene,
    t.cl = if (sex == "Male") sankey_base$male.cl else sankey_base$female.cl,
    n.cl = ncl,
    o.cl = if (sex == "Male") sankey_base$female.cl else sankey_base$male.cl,
    stringsAsFactors = FALSE
  )
  tumor_mat <- make_fit_matrix(trajectory, sex, "Tumor", tumor_features, prediction_ids)
  tumor_df <- tumor_base
  tumor_df$row_split <- factor(paste0(substr(sex, 1, 1), "-T", tumor_df$t.cl),
                               levels = paste0(substr(sex, 1, 1), "-T", 1:4))

  # The published panel computed labels from highlighted pathway/T-TMP genes
  # and then let anno_mark use their actual row positions. Reconstruct that
  # behavior by anchoring labels to the displayed tumor row order directly.
  label_genes <- intersect(mark_genes, tumor_df$gene)
  label_at <- match(label_genes, tumor_df$gene)
  label_order <- order(label_at, seq_along(label_at))
  label_genes <- label_genes[label_order]
  label_at <- label_at[label_order]

  normal_mat0 <- make_fit_matrix(trajectory, sex, "Normal", groups$Gene, prediction_ids)

  # Keep the N-TMP rows in the exact same protein order as the adjacent T-TMP
  # matrix. The normal cluster is an aligned annotation and a secondary ordering
  # variable, not permission to break row-wise tumor-normal correspondence.
  normal_features <- tumor_features
  normal_mat <- normal_mat0[normal_features, , drop = FALSE]
  normal_df <- tumor_df[, c("gene", "n.cl", "t.cl"), drop = FALSE]
  normal_df$row_split <- tumor_df$row_split

  col_df <- data.frame(
    age.class = factor(col_age_class, levels = age_levels),
    stringsAsFactors = FALSE
  )

  rownames(tumor_mat) <- paste0(sex, "_Tumor_", tumor_features)
  rownames(normal_mat) <- paste0(sex, "_Normal_", normal_features)

  list(
    tumor_matrix = tumor_mat,
    normal_matrix = normal_mat,
    tumor_df = tumor_df,
    normal_df = normal_df,
    col_df = col_df,
    label_genes = label_genes,
    label_at = label_at
  )
}

male_panel <- make_sex_panel("Male")
female_panel <- make_sex_panel("Female")

get_green <- function(n) grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[-(1:2)])(n)
cluster_col <- c("1" = "#ef1c1c", "2" = "#3a86c8", "3" = "#2db84d", "4" = "#7a00cc")
age_class_col <- setNames(get_green(4), age_levels)
tissue_col <- c(Tumor = "#FFA500", Normal = "#A020F0")
heatmap_col <- circlize::colorRamp2(c(-1, 0, 1), c("#2B00FF", "white", "#FF2A1A"))

make_tumor_mark_annotation <- function(label_at, labels) {
  keep <- !is.na(label_at) & !is.na(labels)
  label_at <- label_at[keep]
  labels <- labels[keep]
  if (length(label_at) == 0) {
    return(NULL)
  }
  # Labels are placed only on tumor/T-TMP rows. Their vertical anchors are
  # scaled from the full sex-specific tumor matrix order so the reduced
  # normal-overlap display keeps the manuscript-like label heights.
  label_gp <- grid::gpar(
    fontsize = 7,
    col = ifelse(labels %in% c("CNTN1", "MAPT", "L1CAM"), "#377eb8", "black")
  )
  ComplexHeatmap::rowAnnotation(
    labels = ComplexHeatmap::anno_mark(
      at = label_at,
      labels = labels,
      which = "row",
      side = "left",
      labels_gp = label_gp,
      link_width = grid::unit(2, "mm"),
      padding = grid::unit(0.5, "mm"),
      link_gp = grid::gpar(col = "grey30", lwd = 0.7)
    ),
    width = grid::unit(18, "mm"),
    show_annotation_name = FALSE
  )
}

make_group_annotation <- function(values, name) {
  ComplexHeatmap::rowAnnotation(
    group = values,
    col = list(group = cluster_col),
    simple_anno_size = grid::unit(3, "mm"),
    show_annotation_name = FALSE,
    na_col = "gray95"
  )
}

make_top_annotation <- function(type, col_df) {
  ComplexHeatmap::HeatmapAnnotation(
    tissue = factor(rep(type, nrow(col_df)), levels = c("Tumor", "Normal")),
    age.class = col_df$age.class,
    col = list(tissue = tissue_col, age.class = age_class_col),
    annotation_name_gp = grid::gpar(fontsize = 7),
    simple_anno_size = grid::unit(3, "mm"),
    show_legend = TRUE
  )
}

make_tumor_heatmap <- function(panel, sex, show_labels = TRUE) {
  label_annotation <- if (show_labels) make_tumor_mark_annotation(panel$label_at, panel$label_genes) else NULL
  ComplexHeatmap::Heatmap(
    panel$tumor_matrix,
    name = paste0(sex, " Protein TMP"),
    col = heatmap_col,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    use_raster = TRUE,
    raster_quality = 5,
    na_col = "gray95",
    row_gap = grid::unit(5, "mm"),
    column_gap = grid::unit(0.3, "mm"),
    column_split = panel$col_df$age.class,
    row_split = panel$tumor_df$row_split,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = paste0(sex, "\nT-TMP"),
    column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    left_annotation = label_annotation,
    right_annotation = make_group_annotation(panel$tumor_df$t.cl, "T-TMP group"),
    top_annotation = make_top_annotation("Tumor", panel$col_df),
    heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), title = "Protein\nTMP")
  )
}

make_normal_heatmap <- function(panel, sex) {
  ComplexHeatmap::Heatmap(
    panel$normal_matrix,
    name = paste0(sex, " Normal TMP"),
    col = heatmap_col,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    use_raster = TRUE,
    raster_quality = 5,
    na_col = "gray95",
    row_gap = grid::unit(5, "mm"),
    column_gap = grid::unit(0.3, "mm"),
    column_split = panel$col_df$age.class,
    row_split = panel$normal_df$row_split,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = paste0(sex, "\nN-TMP"),
    column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    left_annotation = make_group_annotation(panel$normal_df$n.cl, "N-TMP group"),
    top_annotation = make_top_annotation("Normal", panel$col_df),
    show_heatmap_legend = FALSE
  )
}

message("Drawing Figure 2A with independent male and female row organizations...")
male_ht <- make_tumor_heatmap(male_panel, "Male", show_labels = TRUE) +
  make_normal_heatmap(male_panel, "Male")
female_ht <- make_tumor_heatmap(female_panel, "Female", show_labels = TRUE) +
  make_normal_heatmap(female_panel, "Female")

# Draw the male and female panels separately so the female T-TMP row split is
# controlled by female tumor groups, not by the male panel or normal ordering.
draw_two_panel_heatmap <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 14.5, height = 9.5, useDingbats = FALSE)
  } else {
    grDevices::png(file, width = 14.5, height = 9.5, units = "in", res = 220)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      nrow = 1,
      ncol = 2,
      widths = grid::unit(c(0.5, 0.5), "npc")
    )
  ))

  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  ComplexHeatmap::draw(
    male_ht,
    newpage = FALSE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  )
  grid::upViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  ComplexHeatmap::draw(
    female_ht,
    newpage = FALSE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  )
  grid::upViewport(2)
}

pdf_file <- "figure2a_from_ageTMP.pdf"
png_file <- "figure2a_from_ageTMP.png"
draw_two_panel_heatmap(pdf_file, "pdf")
draw_two_panel_heatmap(png_file, "png")

message("Saved ", pdf_file, " and ", png_file)
