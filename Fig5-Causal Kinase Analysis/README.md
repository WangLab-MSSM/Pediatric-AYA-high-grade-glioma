# Figure 5 — Causal Kinase Analysis

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig5-Causal Kinase Analysis"
```

## Figure Panels

| Analysis/Figures/Tables | Script | Section in the manuscript | Simplified legend |
| :---: | :--- | :--- | :--- |
| Fig 5A | `Fig5_S5.R` | Causal kinase network construction | Topology of a subnetwork illustrating the inferred causal relationship among the 26 kinases with first or second degree causal associations with patient OS in the cDiscovery cohort (ages 0-40). |
| Fig 5B | `Fig5_S5.R` | Causal kinase network construction | Topology of the subnetwork for Module I. |
| Fig 5C | `Fig5_S5.R` | Tumor cell-specific kinase regulation | Heatmap showing estimated tumor cell-specific correlation coefficients among kinase activity scores (KEA3) for the seven genes in Module I. |
| Fig 5D | `Fig5_S5.R` | Tumor cell-specific kinase regulation | Heatmap showing pathways up or down regulated in cancer cells by six Module I kinases as inferred by cell type-specific analysis. |
| Fig 5E | `Fig5_S5.R` | Tumor cell-specific kinase regulation | Volcano plot illustrating tumor cell-specific associations between proteins from the oxidative phosphorylation pathway and CDK8 kinase activities. |
| Fig 5F | `Fig5_S5.R` | Cell line based investigation of precision treatment strategies targeting causal kinases | Heatmap for inferred z scores measuring interaction strength between pairs of kinases in Module I. |
| Fig 5G | `Fig5_S5.R` | Cell line based investigation of precision treatment strategies targeting causal kinases | Diagram illustrating a precision treatment strategy and associated cell line experiments for validation. |
| Fig 5H | `Fig5_S5.R` | Cell line based investigation of precision treatment strategies targeting causal kinases | Bar plots showing concordance between ATM, CDK8, and LCK kinase activity and proliferation scores based on drug inhibition and CRISPR knockout experiments in tumor-derived cell lines. |
| Fig 5I | `Fig5_S5.R` | Cell line based investigation of precision treatment strategies targeting causal kinases | Scatter plots showing the negative correlation between tumor kinase activity scores and cell line drug response z-scores from drug inhibition experiments for ATM, CDK8, and LCK. |
