## ============================================================
## Tumor vs normal protein plots (all cancer groups)
##
## This script generates tumor-versus-normal protein abundance
## plots for selected genes across all tumor samples combined.
##
## Notes:
## - Study-derived clinical variables are harmonized to simple
##   analysis-facing names (e.g., age, Gender) for plotting
##   clarity. Original source columns are retained in the
##   clinical table.
## - Tumor protein data are read directly from the public protein
##   matrix in ../data.
## - Normal reference proteome data are loaded from the ageTMP
##   package reference dataset.
## ============================================================

if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("../ageTMP")) {
  pkgload::load_all("../ageTMP", quiet = TRUE)
} else if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("ageTMP")) {
  pkgload::load_all("ageTMP", quiet = TRUE)
} else if (requireNamespace("ageTMP", quietly = TRUE)) {
  library(ageTMP)
} else {
  stop("Install ageTMP or run this script from the repository root containing ageTMP/.", call. = FALSE)
}

library(ggplot2)
library(ggpubr)
library(dplyr)

##Helper function

mean_collapse_by_gene <- function(data_mat, gene_vec) {
  genes <- unique(gene_vec)
  
  out <- do.call(
    rbind,
    lapply(genes, function(g) {
      apply(
        data_mat[gene_vec == g, , drop = FALSE],
        2,
        mean,
        na.rm = TRUE
      )
    })
  )
  
  rownames(out) <- genes
  out
}
## ------------------------------------------------------------
## Load data
## ------------------------------------------------------------

data_dir <- "../data"

clinical.data <- ageTMP::ageTMP_load_clinical(data_dir)
clinical.data$id <- ageTMP::ageTMP_normalize_sample_ids(clinical.data$id)
if ("sample_id" %in% colnames(clinical.data)) {
  clinical.data$sample_id <- ageTMP::ageTMP_normalize_sample_ids(clinical.data$sample_id)
}

protein.raw <- ageTMP::ageTMP_load_molecular(data_dir, modality = "protein")
protein.split <- ageTMP::ageTMP_split_annotation_matrix(
  protein.raw,
  annotation_cols = 1:4,
  row_id = "ApprovedGeneSymbol"
)
tumor.protein.data <- mean_collapse_by_gene(
  protein.split$matrix,
  protein.split$annotation$ApprovedGeneSymbol
)

normal.reference <- ageTMP::ageTMP_load_normal_reference()
normal.protein.data <- normal.reference$protein$matrix
normal.protein.meta <- normal.reference$protein$sample_metadata

set.seed(123)

## ------------------------------------------------------------
## Harmonize analysis-facing clinical variable names
## ------------------------------------------------------------
## For plotting clarity, create standardized analysis aliases
## from source-specific clinical columns. 

clinical.data$age <- clinical.data$cDisc_age
clinical.data$Gender <- clinical.data$cDisc_Gender

clinical.data$cancer.group <- ifelse(
  clinical.data$cohort == "HOPE",
  as.character(clinical.data$Disc_Cancer_Group),
  "ADULT-GBM"
)

clinical.data$cancer.group <- ifelse(
  !clinical.data$Disc_Cancer_Group %in% c(
    "(DMG) Diffuse Midline Glioma",
    "ADULT-GBM",
    "(HGG) High Grade Glioma (not otherwise specified)"
  ),
  "Other",
  as.character(clinical.data$Disc_Cancer_Group)
)


clinical.data$cancer.group <- dplyr::case_when(
  clinical.data$cancer.group == "(HGG) High Grade Glioma (not otherwise specified)" ~ "HGG-NOS",
  clinical.data$cancer.group == "(DMG) Diffuse Midline Glioma" ~ "DMG",
  clinical.data$cancer.group == "ADULT-GBM" ~ "ADULT-GBM",
  TRUE ~ "Other"
)

## ------------------------------------------------------------
## Restrict to HGG-NOS tumor samples
## ------------------------------------------------------------

nos.ids <- clinical.data[
  clinical.data$cancer.group == "HGG-NOS",
  "id"
]

## ------------------------------------------------------------
## Function: tumor vs normal plot for one gene
## ------------------------------------------------------------

get_tn_plot <- function(mygene) {
  
  ## ---- normal reference samples
  normal.df <- data.frame(
    value = as.numeric(scale(normal.protein.data[match(mygene, rownames(normal.protein.data)), ])),
    type = "normal",
    age = normal.protein.meta[match(colnames(normal.protein.data), normal.protein.meta$ID), ]$Age,
    sex = ifelse(
      normal.protein.meta[match(colnames(normal.protein.data), normal.protein.meta$ID), ]$Gender == "m",
      "Male",
      "Female"
    )
  )
  
  normal.df$age.class <- cut(
    normal.df$age,
    breaks = c(0, 15, 40),
    include.lowest = TRUE,
    labels = c("PED", "AYA")
  )
  normal.df <- normal.df[normal.df$age < 50, ]
  normal.df$value <- as.numeric(scale(normal.df$value))
  
  ## ---- tumor samples (all cancer groups combined)
  tumor.df <- data.frame(
    value = as.numeric(scale(as.numeric(tumor.protein.data[match(mygene, rownames(tumor.protein.data)), ]))),
    type = "tumor",
    age = clinical.data[match(colnames(tumor.protein.data), clinical.data$id), ]$age,
    sex = clinical.data[match(colnames(tumor.protein.data), clinical.data$id), ]$Gender
  )
  
  tumor.df <- tumor.df[tumor.df$age < 50, ]
  tumor.df$age.class <- cut(
    tumor.df$age,
    breaks = c(0, 15, 40),
    include.lowest = TRUE,
    labels = c("PED", "AYA")
  )
  tumor.df$value <- as.numeric(scale(tumor.df$value))
  
  ## ---- combine
  plot.df <- rbind(normal.df, tumor.df)
  plot.df$type <- factor(plot.df$type)
  plot.df$age.class <- factor(plot.df$age.class, levels = c("PED", "AYA"))
  plot.df$sex <- factor(plot.df$sex, levels = c("Male", "Female"))
  plot.df$label <- factor(paste(plot.df$type, plot.df$sex))
  plot.df <- plot.df[!is.na(plot.df$age.class), ]
  
  ## ---- statistical testing
  test.df <- plot.df %>%
    filter(!is.na(value), !is.na(type), !is.na(age.class), !is.na(sex)) %>%
    group_by(sex, age.class) %>%
    filter(n_distinct(type) > 1) %>%
    ungroup()
  
  ## ---- plot
  ggplot(plot.df, aes(x = age.class, y = value, fill = label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(size = 0.5, position = position_jitterdodge(jitter.width = 0.5)) +
    stat_compare_means(
      data = test.df,
      label = "p.signif",
      label.y = 2,
      size = 5,
      fontface = "bold",
      hide.ns = TRUE,
      method = "wilcox.test",
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "")
      )
    ) +
    scale_fill_manual(values = c(
      "normal Female" = "#FFA7A7",
      "normal Male" = "#A3A3FF",
      "tumor Female" = "#CC0303",
      "tumor Male" = "#0707CF"
    )) +
    ggtitle(mygene) +
    facet_grid(cols = vars(sex)) +
    theme_bw() +
    coord_cartesian(ylim = c(-3, 3)) +
    ylab("Protein abundance") +
    xlab("Age class")
}

## ------------------------------------------------------------
## Genes to plot
## ------------------------------------------------------------

mygenes <- c("CNTN1", "MAPT", "L1CAM")
plots <- lapply(mygenes, get_tn_plot)

## ------------------------------------------------------------
## Save figure
## ------------------------------------------------------------

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(output_dir, "Figure2D_protein_tn.pdf"), height = 3, width = 3 * length(mygenes))
print(ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, nrow = 1))
dev.off()
