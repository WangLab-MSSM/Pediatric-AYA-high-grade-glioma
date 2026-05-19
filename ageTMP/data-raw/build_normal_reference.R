## Build normal_reference.rds for ageTMP
##
## The normal/reference data used for AD-TMP trajectory analysis were
## previously stored in the paper-analysis object `for_tadj.RData`.
##
## Source/provenance:
## - DEveLopmental Trajectory Atlas (DELTA), dorsolateral prefrontal cortex
##   (DLPFC)
## - PMID: 30518843
## - URL: http://amp.pharm.mssm.edu/DELTA
##
## This script converts the legacy study object into documented package data.

legacy_file <- Sys.getenv(
  "AGETMP_FOR_TADJ",
  "/Users/lashbn01/Dropbox/HOPE_AYA/for_tadj.RData"
)

if (!file.exists(legacy_file)) {
  stop("Legacy for_tadj.RData file not found: ", legacy_file, call. = FALSE)
}

loaded_objects <- load(legacy_file)

normal_reference <- list(
  protein = list(
    matrix = breen.prot,
    sample_metadata = breen.prot.meta
  ),
  transcript = list(
    matrix = breen.trans,
    sample_metadata = breen.trans.meta
  ),
  brainspan = list(
    annotation = bs.anno,
    matrix = bs.data,
    sample_metadata = bs.meta
  ),
  provenance = list(
    source_name = "DEveLopmental Trajectory Atlas (DELTA), DLPFC",
    pmid = "30518843",
    url = "http://amp.pharm.mssm.edu/DELTA",
    legacy_file = legacy_file,
    legacy_objects = loaded_objects,
    build_date = as.character(Sys.Date())
  )
)

out_dir <- file.path("ageTMP", "inst", "extdata")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

saveRDS(normal_reference, file.path(out_dir, "normal_reference.rds"))
