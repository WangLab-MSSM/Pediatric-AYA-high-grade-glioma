#!/usr/bin/env Rscript

# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai

# Figure 1F final-data analogue. The age classes displayed in the manuscript
# were derived early in the study using then-current data versions. This script
# uses the final repository data tables and temporalCPSA to render the same
# style of AD-TMP heatmap. Protein trajectories reproduce the archived
# visualization input closely, whereas RNA reflects modest source/legacy
# matrix-version differences and phospho reflects a larger final phosphosite
# universe before ApprovedGeneSymbol-level collapse. The source-derived plot is
# therefore intended as a transparent final-data companion analysis, not a
# pixel-identical reproduction of the manuscript panel.
#
# The survival-days track is descriptive. In source-derived mode it defaults to
# the currently available reference clinical table, which has a lower
# observed-survival distribution among event-coded cases than the archived
# figure-preparation annotation; the resulting peak can therefore be lower.
# The smoother span is kept fixed rather than tuned to mimic the manuscript
# visualization.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else ""
script_path <- gsub("~\\+~", " ", script_path)
script_dir <- if (nzchar(script_path)) {
  normalizePath(dirname(script_path), mustWork = TRUE)
} else if (file.exists(file.path(getwd(), "Figure1F_AD_TMP_heatmap_source_derived.R"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(dirname(getwd()), mustWork = TRUE)
}

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
for (modality in c("protein", "rna", "phospho")) {
  env_name <- paste0("AGETMP_", toupper(modality), "_FIT_SAMPLE_REFERENCE")
  if (!nzchar(Sys.getenv(env_name))) {
    env_value <- file.path(script_dir, "legacy_inputs", paste0("figure1f_", modality, "_fit_samples.tsv"))
    do.call(Sys.setenv, stats::setNames(list(env_value), env_name))
  }
}

source(file.path(script_dir, "Figure1F_AD_TMP_heatmap.R"), chdir = FALSE)
