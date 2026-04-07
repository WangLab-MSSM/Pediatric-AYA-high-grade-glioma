# Pediatric AYA High-Grade Glioma (HOPE Study) — Analysis Repository

This repository contains code used to generate all figures and tables for the HOPE pediatric/AYA high-grade glioma study.

---

## Repository structure

Each folder corresponds to a figure in the manuscript and contains the code required to reproduce that analysis.

### Figure-to-folder mapping

| Figures   | Name of the GitHub folder |
|:---------:|:--------------------------|
| Figure 1  | Fig1-Study overview and Landscape of DNA alteration |
| Figure 2  | Fig2-Age-dependent tumor molecular profiles |
| Figure 3  | Fig3-Prognostic markers/pathways based on AD-TMP |
| Figure 4  | Fig4-Causal Kinase Analysis |
| Figure 5  | Fig5-Cis and Trans Regulation |
| Figure 6  | Fig6-Prognostic Subtyping |
| Figure 7  | Fig7-Immune Landscape |

---

## How to use this repository

To reproduce a specific figure:

1. Download the required data (see **Data access** below)
2. Place the data files in the `/data` folder in this repository
3. Open R or RStudio and set your working directory to the repository folder
4. Navigate to the corresponding figure folder
5. Run the scripts in that folder

Each folder contains self-contained analysis code for that figure.

---

## Data access

The analysis-ready datasets are not stored in this repository due to size constraints.

### Data organization (important)

All primary data files are located in the /data directory at the root of this repository.

Scripts are written to load input data from this directory to ensure consistent and reproducible results.

Some figure-specific folders may also contain intermediate or derived data files used for generating individual figures. These are provided for convenience and to facilitate reproduction of specific analyses.

---

### Step 1: Download the data

Download the following files and place them in the `/data` folder (currently empty):

**pediatric_aya_hgg_external_data.rds**  
https://www.dropbox.com/scl/fi/gcd01rm0rzn4iw1e320va/pediatric_aya_hgg_external_data.rds?rlkey=x8obchfrubbdq2qn80aw1rs6w&st=8tn52mln&dl=0

**pediatric_aya_hgg_study_data.rds**  
https://www.dropbox.com/scl/fi/4tvlari5ou4vdq087ma4m/pediatric_aya_hgg_study_data.rds?rlkey=qi8pejja6zyzw766tze5fnubo&st=23nnh1ss&dl=0

If access is restricted, please contact:  
**nicole.tignor@gmail.com**

---

## Requirements

This code was developed using R. Required packages are loaded within each script.