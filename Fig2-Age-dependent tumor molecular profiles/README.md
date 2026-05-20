# Figure 2 — Age-Dependent Tumor Molecular Profiles

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig2-Age-dependent tumor molecular profiles"
```

Scripts that call `library(ageTMP)` or `ageTMP::` require installation of the bundled `ageTMP` package as described in the root `README.md`.

## Figure Panels

| Analysis/Figures/Tables | Script | Section in the manuscript | Simplified legend |
| :---: | :--- | :--- | :--- |
| Figure 2A | `Figure2A_TMP_heatmap.R` | T-TMP Derived Protein Groups Exhibiting Sex Associated Patterns | Heatmap comparing normal and tumor Temporal Molecular Profiles in males and females. |
| Figure 2B | `Figure2B_tn_diff_by_sex.R` | T-TMP Derived Protein Groups Exhibiting Sex Associated Patterns | Summary of pathways enriched in T-TMP protein groups for females and males. |
| Figure 2C | `Figure2C_protein_curveplot.R` | Tumor Normal Differences Revealed by TMPs | Protein T-TMPs and N-TMPs for selected glycoproteins involved in neuron development, with shaded areas representing confidence intervals. |
| Figure 2D | `Figure2D_tn_boxplot.R` | Tumor Normal Differences Revealed by TMPs | Boxplots showing protein abundance distributions in different age groups in normal and tumor data for selected glycoproteins involved in neuron development. |
| Figure 2E | `Figure2E_protein_heatmap.R` | Tumor Normal Differences Revealed by TMPs | Heatmaps displaying associations with abnormality (T/N Distance) and sex bias (M/F Difference) for sex-age groups. |
| Figure 2F | `Figure2F_oxphos_curveplot.R` | Tumor Normal Differences Revealed by TMPs | Proteins in the Oxidative Phosphorylation pathway exhibiting extreme tumor-normal differences. |
| Figure 2G | `Figure2G_neuronal_system_curveplot.R` | Tumor Normal Differences Revealed by TMPs | Proteins in the Neuronal System pathway exhibiting extreme tumor-normal differences. |
| Figure S2E | `FigureS2E_RNA_heatmap.R` | Tumor Normal Differences Revealed by TMPs | RNA heatmap related to tumor-normal age-dependent molecular profiles. |
