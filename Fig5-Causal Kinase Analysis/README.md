# Figure 5 — Causal Kinase Analysis

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig5-Causal Kinase Analysis"
```

Generated PDFs and PNG previews are written to `output/`.

Note: this folder includes a figure-specific `data/` directory containing the
kinase-network support files used by `Figure5_S5_causal_kinase_analysis.R`. Those files are read from the
local Figure 5 folder rather than the shared repository-level `../data/`
directory.

## Figure Panels

| Analysis/Figures/Tables | Script | Output file(s) | Section in the manuscript | Simplified legend |
| :---: | :--- | :--- | :--- | :--- |
| Figure 5A | `Figure5_S5_causal_kinase_analysis.R` | `Figure5A_causal_kinase_network.pdf` | Causal kinase network construction | Topology of a subnetwork illustrating the inferred causal relationship among the 26 kinases with first or second degree causal associations with patient OS in the cDiscovery cohort (ages 0-40). |
| Figure 5B | `Figure5_S5_causal_kinase_analysis.R` | `Figure5B_module_I_causal_kinase_network.pdf` | Causal kinase network construction | Topology of the subnetwork for Module I. |
| Figure 5C | `Figure5_S5_causal_kinase_analysis.R` | `Figure5C_module_I_kinase_correlation_heatmap.pdf` | Tumor cell-specific kinase regulation | Heatmap showing estimated tumor cell-specific correlation coefficients among kinase activity scores (KEA3) for the seven genes in Module I. |
| Figure 5D | `Figure5_S5_causal_kinase_analysis.R` | No standalone output from this script | Tumor cell-specific kinase regulation | Heatmap showing pathways up or down regulated in cancer cells by six Module I kinases as inferred by cell type-specific analysis. |
| Figure 5E | `Figure5_S5_causal_kinase_analysis.R` | `Figure5E_CDK8_oxphos_volcano.pdf` | Tumor cell-specific kinase regulation | Volcano plot illustrating tumor cell-specific associations between proteins from the oxidative phosphorylation pathway and CDK8 kinase activities. |
| Figure 5F | `Figure5_S5_causal_kinase_analysis.R` | `Figure5F_module_I_kinase_interaction_heatmap.pdf` | Cell line based investigation of precision treatment strategies targeting causal kinases | Heatmap for inferred z scores measuring interaction strength between pairs of kinases in Module I. |
| Figure 5G | `Figure5_S5_causal_kinase_analysis.R` | No standalone output from this script | Cell line based investigation of precision treatment strategies targeting causal kinases | Diagram illustrating a precision treatment strategy and associated cell line experiments for validation. |
| Figure 5H | `Figure5_S5_causal_kinase_analysis.R` | `Figure5H_ATM_CDK8_LCK_validation_barplot.pdf` | Cell line based investigation of precision treatment strategies targeting causal kinases | Bar plots showing concordance between ATM, CDK8, and LCK kinase activity and proliferation scores based on drug inhibition and CRISPR knockout experiments in tumor-derived cell lines. |
| Figure 5I | `Figure5_S5_causal_kinase_analysis.R` | `Figure5I_ATM_CDK8_LCK_drug_response_scatter.pdf` | Cell line based investigation of precision treatment strategies targeting causal kinases | Scatter plots showing the negative correlation between tumor kinase activity scores and cell line drug response z-scores from drug inhibition experiments for ATM, CDK8, and LCK. |

## Supplementary Panels

| Analysis/Figures/Tables | Script | Output file(s) | Simplified legend |
| :---: | :--- | :--- | :--- |
| Figure S5B | `Figure5_S5_causal_kinase_analysis.R` | `FigureS5B_module_II_causal_kinase_network.pdf` | Topology of the subnetwork for Module II. |
| Figure S5C | `Figure5_S5_causal_kinase_analysis.R` | `FigureS5C_male_module_I_causal_kinase_network.pdf` | Male-stratified topology of the Module I subnetwork. |
| Figure S5D | `Figure5_S5_causal_kinase_analysis.R` | `FigureS5D_male_module_II_causal_kinase_network.pdf` | Male-stratified topology of the Module II subnetwork. |
| Figure S5E | `Figure5_S5_causal_kinase_analysis.R` | `FigureS5E_module_II_kinase_correlation_heatmap.pdf` | Tumor cell-specific correlation heatmap for Module II kinase activity scores. |
| Figure S5G | `Figure5_S5_causal_kinase_analysis.R` | `FigureS5G_CDK8_survival_curve.pdf` | Survival curve based on CDK8 expression in the reference cohort. |
| Figure S5I | `Figure5_S5_causal_kinase_analysis.R` | `FigureS5I_FLT1_MAP2K1_RAF1_validation_barplot.pdf` | Concordance between FLT1, MAP2K1, and RAF1 kinase activity and proliferation scores. |
| Figure S5J | `Figure5_S5_causal_kinase_analysis.R` | `FigureS5J_FLT1_MAP2K1_RAF1_drug_response_scatter.pdf` | Drug response scatter plots for FLT1, MAP2K1, and RAF1. |

## Rendered Output Checklist

Running `Figure5_S5_causal_kinase_analysis.R` writes panel-labeled PDFs to
`output/`. PNG files with the same basename may be created for visual
inspection.

| Panel label | PDF output | Optional preview |
| :---: | :--- | :--- |
| Figure 5A | `Figure5A_causal_kinase_network.pdf` | `Figure5A_causal_kinase_network.png` |
| Figure 5B | `Figure5B_module_I_causal_kinase_network.pdf` | `Figure5B_module_I_causal_kinase_network.png` |
| Figure 5C | `Figure5C_module_I_kinase_correlation_heatmap.pdf` | `Figure5C_module_I_kinase_correlation_heatmap.png` |
| Figure 5E | `Figure5E_CDK8_oxphos_volcano.pdf` | `Figure5E_CDK8_oxphos_volcano.png` |
| Figure 5F | `Figure5F_module_I_kinase_interaction_heatmap.pdf` | `Figure5F_module_I_kinase_interaction_heatmap.png` |
| Figure 5H | `Figure5H_ATM_CDK8_LCK_validation_barplot.pdf` | `Figure5H_ATM_CDK8_LCK_validation_barplot.png` |
| Figure 5I | `Figure5I_ATM_CDK8_LCK_drug_response_scatter.pdf` | `Figure5I_ATM_CDK8_LCK_drug_response_scatter.png` |
| Figure S5B | `FigureS5B_module_II_causal_kinase_network.pdf` | `FigureS5B_module_II_causal_kinase_network.png` |
| Figure S5C | `FigureS5C_male_module_I_causal_kinase_network.pdf` | `FigureS5C_male_module_I_causal_kinase_network.png` |
| Figure S5D | `FigureS5D_male_module_II_causal_kinase_network.pdf` | `FigureS5D_male_module_II_causal_kinase_network.png` |
| Figure S5E | `FigureS5E_module_II_kinase_correlation_heatmap.pdf` | `FigureS5E_module_II_kinase_correlation_heatmap.png` |
| Figure S5G | `FigureS5G_CDK8_survival_curve.pdf` | `FigureS5G_CDK8_survival_curve.png` |
| Figure S5I | `FigureS5I_FLT1_MAP2K1_RAF1_validation_barplot.pdf` | `FigureS5I_FLT1_MAP2K1_RAF1_validation_barplot.png` |
| Figure S5J | `FigureS5J_FLT1_MAP2K1_RAF1_drug_response_scatter.pdf` | `FigureS5J_FLT1_MAP2K1_RAF1_drug_response_scatter.png` |

Figure 5D and Figure 5G are listed in the manuscript panel map above, but they
do not currently have standalone plot outputs in this script.
