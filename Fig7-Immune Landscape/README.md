# Figure 7 — Immune Landscape

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig7-Immune Landscape"
```

Assembled Figure 7 PDFs and PNG previews are stored in `output/`.

## Figure Panels

| Analysis/Figures/Tables | Script | Output file(s) | Section in the manuscript | Simplified legend |
| :---: | :--- | :--- | :--- | :--- |
| Figure 7A | `Figure7_BayesDeBulk_run.R` provides the BayesDeBulk setup/run; assembled figure output provided separately | Included in assembled Figure 7 output | Immune signaling in PED/AYA HGG tumors | Heatmap showing estimated cell type fractions (z-scored) by BayesDeBulk for the cDiscovery tumors. |
| Figure 7B | `Figure7_BayesDeBulk_run.R` provides the BayesDeBulk setup/run; assembled figure output provided separately | Included in assembled Figure 7 output | Overall survival outcome differences across Protein Subtypes | Boxplots showing estimated fractions of T cells and SEPP1 Mo TAM stratified by age classes. |
| Figure 7C | `Figure7_BayesDeBulk_run.R` provides the BayesDeBulk setup/run; assembled figure output provided separately | `MainFigure7_CD.pdf` | Immune signaling in PED/AYA HGG tumors | Boxplots showing estimated fractions of T cells and SEPP1 Mo TAM stratified by protein clusters. |
| Figure 7D | `Figure7_BayesDeBulk_run.R` provides the BayesDeBulk setup/run; assembled figure output provided separately | `MainFigure7_CD.pdf` | Immune signaling in PED/AYA HGG tumors | Survival curves comparing female patients (ages 0-40) with low SEPP1 Mo TAM and high T-cell infiltration to those with high SEPP1 Mo TAM and low T-cell infiltration in the cDiscovery cohort. |
| Figure 7E | `Figure7_BayesDeBulk_run.R` provides the BayesDeBulk setup/run; assembled figure output provided separately | Included in assembled Figure 7 output | Immune signaling in PED/AYA HGG tumors | Survival curves from the sex-adjusted Cox regression model comparing patients with low SEPP1 Mo TAM and high T-cell infiltration to those with high SEPP1 Mo TAM and low T-cell infiltration in the Validation cohort. |
| Figure 7F | `Figure7_BayesDeBulk_run.R` provides the BayesDeBulk setup/run; assembled figure output provided separately | Included in assembled Figure 7 output | Immune signaling in PED/AYA HGG tumors | Boxplots showing estimated fractions of T cells and SEPP1 Mo TAM stratified by NF1 mutation status within the cDiscovery cohort. |
| Figure 7G | `Figure7_BayesDeBulk_run.R` provides the BayesDeBulk setup/run; assembled figure output provided separately | Included in assembled Figure 7 output | Immune signaling in PED/AYA HGG tumors | Boxplots showing estimated fractions of T cells and SEPP1 Mo TAM stratified by H3-3A mutation status within the cDiscovery cohort. |
| Figure 7H | `Figure7_BayesDeBulk_run.R` provides the BayesDeBulk setup/run; assembled figure output provided separately | Included in assembled Figure 7 output | Immune signaling in PED/AYA HGG tumors | Heatmap showing associations between estimated T-cell and SEPP1 Mo TAM fractions and immune checkpoint gene/protein abundances across sex and age groups. |

## Current Script Status

`Figure7_BayesDeBulk_run.R` performs the BayesDeBulk preprocessing and model run
using `Data.rda` and `Gene_signature.rda`, producing `pi.post` in memory. This
script represents the computational setup for the Figure 7 immune deconvolution
analysis.

The current repository also includes assembled Figure 7 C/D outputs:

- `MainFigure7_CD.pdf`
- `SFigure7_CD.pdf`

These assembled PDFs collect the available Figure 7 C/D panels rather than
writing each panel as a separate file from the BayesDeBulk run script. If
panel-wise Figure 7 plotting code is added later, output PDFs should be written
to `output/` and should include the panel label in the filename.
