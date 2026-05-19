# Figure 4 — Prognostic Markers and Pathways Based on AD-TMP

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig4-Prognostic markers:pathways based on AD-TMP"
```

Scripts that call `library(ageTMP)` or `ageTMP::` require installation of the bundled `ageTMP` package as described in the root `README.md`.

## Figure Panels

| Analysis/Figures/Tables | Script | Section in the manuscript | Simplified legend |
| :--- | :--- | :--- | :--- |
| Fig 4A | `figure4_sex_specific_cox_models.R` | Validating a prognostic score based on wild-type IDH1 and IDH2 proteins | Survival plots comparing OS of male, H3-wild-type patients with high or low IDH1/2 scores in the cDiscovery cohort using TMT or SRM protein data. |
| Fig 4B | `figure4b_survival_assoc_glyco.R` | Glycosylation associated with overall survival outcomes | Heatmap illustrating OS-associated glycopeptides from OS-associated pathways. |
| Fig 4C | `figure4c_survival_landscape.R` | Age-related neurodevelopment and tumor progression via CPSA | Heatmap illustrating overlap among genes with significant OS associations across sex and age groups from CPSA using global protein or RNA data. |
| Fig 4D | `figure4d_survival_assoc_heatmap.R` | OS-associated pathways in different sex-age strata | Heatmaps illustrating concordant OS association results from CPSA using global protein data. |
| Fig 4E; Fig S4D; Fig S4E; Fig S4F; Fig 4I; Fig S4I; STable4 | `figure4_sex_specific_cox_models.R` | Construction and validation of protein-based OS risk scores | Fits sex-specific ridge-penalized Cox models in the cDiscovery cohort, derives frozen male and female prognostic score coefficients, generates coefficient plots, model-adjusted survival plots, HR tables, and downstream-ready objects. |
| Fig 4G; Fig 4H; Fig S4 combined-score panels | `figure4_combined_frozen_score_models.R` | Construction and validation of protein-based OS risk scores | Applies frozen male- and female-derived prognostic score coefficients to cDiscovery and Validation cohorts, fits combined Cox models containing both scores, and generates multivariable HR plots/tables plus adjusted survival plots. |
| Fig 4J and related validation outputs | `figure4_validation_combined_score_model.R` | Construction and validation of protein-based OS risk scores | Validates selected prognostic marker and combined-score models in the independent Validation cohort. |

## Run Order Notes

`figure4_combined_frozen_score_models.R` uses outputs created by `figure4_sex_specific_cox_models.R`, so run `figure4_sex_specific_cox_models.R` first when reproducing the frozen-score analyses.
