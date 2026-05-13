## ============================================================
## Build study data object for Pediatric AYA HGG
##
## This script builds the study-derived analysis object from:
##   1) DCC molecular data files
##   2) manuscript Supplementary Tables
##
## Output:
##   data/pediatric_aya_hgg_study_data.rds
## ============================================================

base_dir <- "data"

read_tsv_file <- function(fname, dir = base_dir) {
  path <- file.path(dir, fname)
  message("Reading: ", fname)
  out <- read.delim(path, sep = "\t", check.names = FALSE)
  message("  dim = ", paste(dim(out), collapse = " x "))
  out
}

split_anno_data <- function(df, anno_cols, row_id = NULL) {
  anno <- df[, anno_cols, drop = FALSE]
  data <- df[, -anno_cols, drop = FALSE]
  
  if (!is.null(row_id) && row_id %in% colnames(anno)) {
    rownames(data) <- anno[[row_id]]
  }
  
  list(anno = anno, data = data)
}


## ------------------------------------------------------------
## Read DCC source files
## ------------------------------------------------------------




file_names= c(
  "cDisc_clinical_data_04032026.tsv",
  "cDisc_gene_location_10232023.tsv",
  "cDisc_CNV_coding_10252023.tsv",
  "cDisc_fusion_feature_10192023.tsv",
  "cDisc_fusion_gene_10192023.tsv",
  "cDisc_methylation_gene_10192023.tsv",
  "cDisc_mutation_10192023.tsv",
  "cDisc_phosphosite_imputed_data_ischemia_removed_motif_11032023.tsv",
  "cDisc_proteome_imputed_data_09152023.tsv",
  "cDisc_rna_coding_10192023.tsv",
  "Disc_full_mutation_data_100224.tsv",
  "Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv"
 
)

file_list <- lapply(file_names, read_tsv_file)
names(file_list) <- file_names

## ------------------------------------------------------------
## Extract study datasets from DCC
## ------------------------------------------------------------

gene_annotation <- file_list[["cDisc_gene_location_10232023.tsv"]]

glyco <- split_anno_data(
  file_list[["Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv"]],
  anno_cols = 1:4
)

cnv <- split_anno_data(
  file_list[["cDisc_rna_coding_10192023.tsv"]],
  anno_cols = 1
)


fusion_feature <- split_anno_data(
  file_list[["cDisc_fusion_feature_10192023.tsv"]],
  anno_cols = 1
)

fusion_gene <- split_anno_data(
  file_list[["cDisc_fusion_gene_10192023.tsv"]],
  anno_cols = 1
)

phospho <- split_anno_data(
  file_list[["cDisc_phosphosite_imputed_data_ischemia_removed_motif_11032023.tsv"]],
  anno_cols = 1:9
)

protein <- split_anno_data(
  file_list[["cDisc_proteome_imputed_data_09152023.tsv"]],
  anno_cols = 1:4
)

rna <- split_anno_data(
  file_list[["cDisc_rna_coding_10192023.tsv"]],
  anno_cols = 1:3,
  row_id = "ApprovedGeneSymbol"
)

mutation <- split_anno_data(
  file_list[["cDisc_mutation_10192023.tsv"]],
  anno_cols = 1,
  row_id = "ApprovedGeneSymbol"
)



## ------------------------------------------------------------
## Load study-derived supplementary data
## ------------------------------------------------------------

clinical_data <- read.delim(
  "data/cDisc_clinical_data_04032026.tsv",
 sep="\t"
)

clinical_data$sample_id=ifelse(grepl("^C",clinical_data$id),gsub("-",".",clinical_data$id,),paste0("X",gsub("-",".",clinical_data$id)))



## ------------------------------------------------------------
## Build study data object
## ------------------------------------------------------------

study_data <- list(
  metadata = list(
    project = "Pediatric AYA High-Grade Glioma",
    build_date = as.character(Sys.Date())
  ),
  
  clinical = clinical_data,
  
  
  gene_annotation = gene_annotation,
  
  glyco = glyco,
  cnv = cnv,
  fusion_feature = fusion_feature,
  fusion_gene = fusion_gene,
  
  rna = rna,
  protein = protein,
  
  phospho = phospho,
 
  
  mutation = mutation
)


## ------------------------------------------------------------
## Save
## ------------------------------------------------------------

saveRDS(
  study_data,
  file = "data/pediatric_aya_hgg_study_data.rds"
)


##


