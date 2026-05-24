#!/usr/bin/env Rscript

# Figure 1F final-data analogue. The age classes displayed in the manuscript
# were derived early in the study using then-current data versions. This script
# uses the final repository data tables and temporalCPSA to render the same
# style of AD-TMP heatmap. In this source-derived version, temporalCPSA
# functions generate the AD-TMP/clustme-style matrix from source data before
# plotting, illustrating the age-class distinctions in the final dataset. It is
# intended as a transparent companion analysis, not a pixel-identical
# reproduction of the manuscript panel.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else getwd()
script_path <- gsub("~\\+~", " ", script_path)
script_dir <- normalizePath(dirname(script_path), mustWork = TRUE)

if (!nzchar(Sys.getenv("AGETMP_FIGURE1F_MODE"))) {
  Sys.setenv(AGETMP_FIGURE1F_MODE = "source")
}
if (!nzchar(Sys.getenv("AGETMP_OUTPUT_DIR"))) {
  Sys.setenv(
    AGETMP_OUTPUT_DIR = file.path(
      script_dir,
      "output",
      "figure1f_source_derived_legacy_samples"
    )
  )
}
if (!nzchar(Sys.getenv("AGETMP_SOURCE_SAMPLE_REFERENCE_MATRIX"))) {
  Sys.setenv(
    AGETMP_SOURCE_SAMPLE_REFERENCE_MATRIX = file.path(
      script_dir,
      "legacy_inputs",
      "figure1f_ad_tmp_legacy_clustme.tsv"
    )
  )
}

source(file.path(script_dir, "Figure1F_AD_TMP_heatmap.R"), chdir = FALSE)
