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
## - Tumor protein data are taken from the study data object.
## - Normal reference proteome data are taken from the reference
##   data object.
## ============================================================

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

study.data <- readRDS("pediatric_aya_hgg_study_data.rds")
external.data <- readRDS("pediatric_aya_hgg_external_data.rds")

clinical.data <- study.data$clinical



tumor.protein.data <- mean_collapse_by_gene(study.data$protein$data,
                                            study.data$protein$anno$ApprovedGeneSymbol)


normal.protein.data <- external.data$normal_brain_reference$breen.prot
normal.protein.meta <- external.data$normal_brain_reference$breen.prot.meta

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

pdf("Figure2D_protein_tn.pdf", height = 3, width = 3 * length(mygenes))
ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, nrow = 1)
dev.off()