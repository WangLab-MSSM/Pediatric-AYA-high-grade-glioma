# Figure 6 — Prognostic Subtyping

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig6-Prognostic Subtyping"
```

## Figure Panels

| Analysis/Figures/Tables | Script | Section in the manuscript | Simplified legend |
| :---: | :--- | :--- | :--- |
| Fig 6A | `figure6d_s6b_cluster_annotation_corrplots.R` | Protein Pathway Activities Reveal Distinct Subtypes of High Grade Glioma | Heatmap showing distinct protein pathway activities across protein clusters identified in the Discovery cohort (ages 0-62). |
| Fig 6B | `figure6d_s6b_cluster_annotation_corrplots.R` | Protein Pathway Activities Reveal Distinct Subtypes of High Grade Glioma | Heatmap illustrating protein pathway ssGSEA scores for the same pathways as shown in A within the validation data. |
| Fig 6C | `figure6d_s6b_cluster_annotation_corrplots.R` | Protein Pathway Activities Reveal Distinct Subtypes of High Grade Glioma | Bubble heatmap illustrating correlations between average pathway scores of each protein cluster in the cDiscovery cohort versus the Validation cohort. |
| Fig 6D | `figure6d_s6b_cluster_annotation_corrplots.R` | Protein Pathway Activities Reveal Distinct Subtypes of High Grade Glioma | Bubble heatmap illustrating distributions of H3-3A mutated tumors across protein clusters. |
| Fig 6E | `figure6m_s6o_s6e_s6j_6i_cd276_itgav_pla2g4a_subtype_analysis.R` | Causal kinases and CAR-T markers associated with Protein Subtypes | Boxplots of ATM, LCK, NTRK2, and CDK8 kinase activities in the cDiscovery cohort (ages 0-40), stratified by protein clusters. |
| Fig 6F | `figure6m_s6o_s6e_s6j_6i_cd276_itgav_pla2g4a_subtype_analysis.R` | Causal kinases and CAR-T markers associated with Protein Subtypes | Boxplot of CD276 (B7-H3) protein abundance in the cDiscovery cohort (ages 0-40), stratified by sex and protein clusters. |
| Fig 6G | `Figure6G_subtype_survival_analysis.R` | Overall survival outcome differences across Protein Subtypes | Survival curves comparing F2, M2, and other clusters within the cDiscovery cohort (ages 0-40). |
| Fig 6H | `figure6m_s6o_s6e_s6j_6i_cd276_itgav_pla2g4a_subtype_analysis.R` | Overall survival outcome differences across Protein Subtypes | Volcano plot comparing protein abundances between F2 and M2 for genes in the Interferon Gamma Response pathway. |
| Fig 6I | `figure6m_s6o_s6e_s6j_6i_cd276_itgav_pla2g4a_subtype_analysis.R` | Overall survival outcome differences across Protein Subtypes | Boxplots illustrating protein and phosphosite abundance of PLA2G4A across sexes and protein clusters. |
| Fig 6J | `figure6m_s6o_s6e_s6j_6i_cd276_itgav_pla2g4a_subtype_analysis.R` | Overall survival outcome differences across Protein Subtypes | Scatterplot showing correlations between LCK kinase activity and PLA2G4A protein expression within each sex in the cDiscovery cohort (ages 0-40). |
| Fig 6K | `figure6k_sex_biased_glycopeptide_pathway_enrichment.R` | Overall survival outcome differences across Protein Subtypes | Bar plot illustrating pathway enrichment results for comparing male and female tumors in C2 based on glycopeptide abundance data. |
| Fig 6L | `figure6k_sex_biased_glycopeptide_pathway_enrichment.R` | Overall survival outcome differences across Protein Subtypes | Volcano plot comparing glycopeptide abundances between F2 and M2 for genes in the Epithelial Mesenchymal Transition pathway. |
| Fig 6M | `figure6m_s6o_s6e_s6j_6i_cd276_itgav_pla2g4a_subtype_analysis.R` | Overall survival outcome differences across Protein Subtypes | Boxplot illustrating abundance of the ITGAV glycopeptide with an N-linked sialylated glycan (N3H5F1S1G0) across sex and protein clusters. |
| Fig 6N | `figure6n_itgav_sialylated_glycopeptide_pdl1_correlation.R` | Overall survival outcome differences across Protein Subtypes | Scatterplot showing correlations between global protein abundance of PD-L1 and abundance of the ITGAV glycopeptide N3H5F1S1G0 in female and male tumors. |
