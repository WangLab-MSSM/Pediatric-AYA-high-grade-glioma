# Pediatric AYA High-Grade Glioma (HOPE Study) — Analysis Repository

This repository contains code used to generate all figures and tables for the HOPE pediatric/AYA high-grade glioma study.

-------------------------------------------------------------------------------

REPRODUCIBILITY OVERVIEW

1. Download supplementary tables from the journal website
2. Download raw data from PDC
3. Build the study data object
4. Consult the README for each figure of interest

-------------------------------------------------------------------------------

REPOSITORY STRUCTURE

Each folder corresponds to a figure in the manuscript and contains the code
required to reproduce that analysis.

Figure-to-folder mapping:

Figure 1
  Fig1-Study overview and Landscape of DNA alteration

Figure 2
  Fig2-Age-dependent tumor molecular profiles

Figure 3
  Fig3-Cis and Trans Regulation

Figure 4
  Fig4-Prognostic markers/pathways based on AD-TMP

Figure 5
  Fig5-Causal Kinase Analysis

Figure 6
  Fig6-Prognostic Subtyping

Figure 7
  Fig7-Immune Landscape

-------------------------------------------------------------------------------

DATA AVAILABILITY AND SETUP

STEP 1: Download supplementary tables from the journal website

Download Supplementary Tables STable1–STable7 associated with the manuscript.

7 supplementary table files should be downloaded and placed in:

data/

Do not create subdirectories.
Do not rename files.

-------------------------------------------------------------------------------

STEP 2: Download raw data from PDC

1. Go to:
   https://proteomic.datacommons.cancer.gov/pdc/

2. Search for:
   PDC000497
   PDC000498
   PDC000499

3. Download all associated data tables.

14 data files should be downloaded and placed in:

data/

Do not create subdirectories.
Do not rename files.

-------------------------------------------------------------------------------
STEP 3 (OPTIONAL): Build master study data object

Some legacy scripts use a precompiled study data object for convenience.

Most updated figure scripts now read directly from files in the data/
directory and do not require this object.

To generate the object for compatibility with older scripts, run:

source("build_study_data00.R")

This will generate:

data/pediatric_aya_hgg_analysis_data.rds

See individual figure folder READMEs for script-specific requirements.

-------------------------------------------------------------------------------

STEP 4: Run figure-specific analyses

After generating the study data object:

1. Navigate to the folder corresponding to the figure of interest.
2. Follow the instructions in that folder’s README.md.
3. Run the scripts to reproduce the analysis.

Each figure folder contains self-contained code and any additional instructions
required.

-------------------------------------------------------------------------------

REQUIREMENTS

This code was developed using R.
Required packages are loaded within each script.