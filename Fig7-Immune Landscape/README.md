# Figure 7 — Immune Landscape

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig7-Immune Landscape"
```

## Figure Panels

| Analysis/Figures/Tables | Script | Section in the manuscript | Simplified legend |
| :---: | :--- | :--- | :--- |
| Fig 7A | `BayesDeBulk_run.R` | Immune signaling in PED/AYA HGG tumors | Heatmap showing estimated cell type fractions (z-scored) by BayesDeBulk for the cDiscovery tumors. |
| Fig 7B | `BayesDeBulk_run.R` | Overall survival outcome differences across Protein Subtypes | Boxplots showing estimated fractions of T cells and SEPP1 Mo TAM stratified by age classes. |
| Fig 7C | `BayesDeBulk_run.R` | Immune signaling in PED/AYA HGG tumors | Boxplots showing estimated fractions of T cells and SEPP1 Mo TAM stratified by protein clusters. |
| Fig 7D | `BayesDeBulk_run.R` | Immune signaling in PED/AYA HGG tumors | Survival curves comparing female patients (ages 0-40) with low SEPP1 Mo TAM and high T-cell infiltration to those with high SEPP1 Mo TAM and low T-cell infiltration in the cDiscovery cohort. |
| Fig 7E | `BayesDeBulk_run.R` | Immune signaling in PED/AYA HGG tumors | Survival curves from the sex-adjusted Cox regression model comparing patients with low SEPP1 Mo TAM and high T-cell infiltration to those with high SEPP1 Mo TAM and low T-cell infiltration in the Validation cohort. |
| Fig 7F | `BayesDeBulk_run.R` | Immune signaling in PED/AYA HGG tumors | Boxplots showing estimated fractions of T cells and SEPP1 Mo TAM stratified by NF1 mutation status within the cDiscovery cohort. |
| Fig 7G | `BayesDeBulk_run.R` | Immune signaling in PED/AYA HGG tumors | Boxplots showing estimated fractions of T cells and SEPP1 Mo TAM stratified by H3-3A mutation status within the cDiscovery cohort. |
| Fig 7H | `BayesDeBulk_run.R` | Immune signaling in PED/AYA HGG tumors | Heatmap showing associations between estimated T-cell and SEPP1 Mo TAM fractions and immune checkpoint gene/protein abundances across sex and age groups. |
