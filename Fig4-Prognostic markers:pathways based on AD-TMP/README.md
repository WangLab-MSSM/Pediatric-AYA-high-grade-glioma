### Code for Figure 4 analyses, figures, and tables

This directory contains scripts used to generate Figure 4, selected Figure S4 panels, and related supplementary tables. Some scripts generate multiple outputs, including main figure panels, supplementary figure panels, hazard-ratio tables, model coefficient tables, and intermediate model objects.

### Summary of Figure 4 analysis scripts and outputs

| Figures / Tables | Script | Section in the manuscript | Simplified legend |
| :--- | :--- | :--- | :--- |
| Fig 4A | `figure4_idh1_idh2_score_survival.R` | Validating a prognostic score based on wild-type IDH1 and IDH2 proteins | Survival plots comparing OS of male, H3-wild-type patients with high or low IDH1/2 scores in the cDiscovery cohort using either TMT or SRM protein data. |
| Fig 4B | `figure4_glycopeptide_os_heatmap.R` | Glycosylation associated with overall survival outcomes | Heatmap illustrating OS-associated glycopeptides from OS-associated pathways. |
| Fig 4C | `figure4_cpsa_gene_overlap_heatmap.R` | Age-related neurodevelopment and tumor progression via CPSA | Heatmap illustrating overlap among genes with significant OS associations across sex and age groups from CPSA using global protein or RNA data. |
| Fig 4D | `figure4_cpsa_pathway_and_protein_heatmaps.R` | OS-associated pathways in different sex-age strata | Heatmaps illustrating concordant OS association results from CPSA using global protein data. The left heatmap shows pathway enrichment results; the right heatmap shows selected protein associations. |
| Fig 4E; Fig S4D; Fig S4E; Fig S4F; Fig 4I; Fig S4I; STable4 | `figure4_sex_specific_cox_models.R` | Construction and validation of protein-based OS risk scores | Fits sex-specific ridge-penalized Cox models in the cDiscovery cohort, derives frozen male and female prognostic score coefficients, generates coefficient plots, model-adjusted survival plots, HR tables, and saves downstream-ready objects. |
| Fig 4G; Fig 4H; Fig S4 combined-score panels | `figure4_combined_frozen_score_models.R` | Construction and validation of protein-based OS risk scores | Applies frozen male- and female-derived prognostic score coefficients to cDiscovery and Validation cohorts, fits combined Cox models containing both scores, and generates multivariable HR plots/tables plus adjusted survival plots. |
| Fig 4J and related validation outputs | `figure4_validation_MORC2_PHKG2.R` | Construction and validation of protein-based OS risk scores | Validates male-derived MORC2 and PHKG2 prognostic markers in the independent Validation cohort using continuous Cox models and median-dichotomized adjusted survival plots. |