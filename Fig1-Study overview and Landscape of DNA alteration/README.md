# Figure 1 — Study Overview and Landscape of DNA Alteration

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig1-Study overview and Landscape of DNA alteration"
```

Generated PDFs, PNG previews, and tables are written to `output/`.

## Figure Panels

| Analysis/Figures/Tables | Script | Section in the manuscript | Simplified legend |
| :---: | :--- | :--- | :--- |
| Figure 1A | `Figure1A_clinical_overview.R` | Proteogenomics profiling of the Discovery cohort | Summary of the demographic and clinical annotations of the Discovery study cohort (N=93). |
| Figure 1B | `Figure1B_barplot.R` | Methylation and Genomically informed Cancer Subtype Groups | Distribution of cancer groups stratified by PED and AYA. |
| Figure 1C | `Figure1C_oncplot.R` | Age dependent Mutational and CNV Landscape of HGG | Oncoprint plot for the frequent (>5%) mutations detected in the Discovery Cohort (N=89). |
| Figure 1D | `Figure1E_reference_assoc.R` | Age dependent Mutational and CNV Landscape of HGG | Mutation frequency profiles along age for two groups of genes with distinct mutation patterns based on the Reference cohort data. |
| Figure 1E | `Figure1E_reference_assoc.R` | Age dependent Mutational and CNV Landscape of HGG | Association of overall survival with somatic mutation events across age groups in the Reference cohort, evaluated using multivariate Cox regression models. |
| Figure 1F | `Figure1F_AD_TMP_heatmap.R`; `Figure1F_AD_TMP_heatmap_source_derived.R` | Transcriptomic and Proteomic Profiles Defining Adolescent versus Young Adult Age Groups | The first script reproduces the manuscript panel using the archived manuscript visualization input; the second renders a source-derived analogue from the final repository data tables using `temporalCPSA`. |
| CNV landscape | `Figure1_CNV_landscape.R` | Age dependent Mutational and CNV Landscape of HGG | Copy-number landscape overview for the Discovery cohort. |
| Figure 1G | `Figure1G_reference_survival_assoc.R` | Transcriptomic and Proteomic Profiles Defining Adolescent versus Young Adult Age Groups | Kaplan-Meier overall survival curves for different age groups in the Reference cohort. |
| Figure S1A | `FigureS1A_clinical_overview.R` | Supplementary cohort annotations | Related cohort annotation summary. |

## Run Order

Most Figure 1 scripts can be run independently.

### Figure 1F

Figure 1F has two companion scripts. The manuscript reproduction script renders
the published panel using the archived visualization input generated during
figure preparation. Both Figure 1F scripts require `temporalCPSA`.

```bash
Rscript Figure1F_AD_TMP_heatmap.R
```

The age-class structure shown in Figure 1F was derived early in the study using
the data versions available at that stage of manuscript development. To show the
same visualization concept on the final repository data tables, run the
source-derived analogue:

```bash
Rscript Figure1F_AD_TMP_heatmap_source_derived.R
```

The source-derived analogue uses `temporalCPSA` functions to generate the
AD-TMP matrix from the final source tables before plotting. It is
intended to illustrate the age-class distinctions in the final dataset and is
not expected to be pixel-identical to the manuscript panel.

The script reads manuscript data from `../data/` and writes generated files to
`output/`.
