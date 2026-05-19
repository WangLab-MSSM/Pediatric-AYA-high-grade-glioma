# ageTMP

`ageTMP` is a reproducibility-first R package for modeling and interpreting age-associated molecular trajectories in tumors relative to developmental or reference populations.

The package scope is broader than Cross-Population Survival Analysis alone. CPSA is treated as one downstream analytical module within the broader `ageTMP` framework.

## Current Focus

The first milestone is direct manuscript reproducibility:

```r
install.packages("remotes")
remotes::install_github("lashbroz/ageTMP")

# Example manuscript workflow pieces
ageTMP::ageTMP_data_sources("data")
ageTMP::ageTMP_load_normal_reference()
```

The package is being organized around these modules:

- public data loading and sample harmonization from manuscript source tables;
- age-dependent tumor molecular phenotype and trajectory generation;
- tumor versus normal/reference trajectory comparison, including Figure 2-style analyses;
- trajectory divergence summaries and visualization;
- clustering of age-dependent trajectory patterns;
- downstream survival modeling, including CPSA.

Normal/reference DLPFC developmental data used for trajectory analyses are documented package extdata derived from the DEveLopmental Trajectory Atlas (DELTA), PMID: 30518843, available at <http://amp.pharm.mssm.edu/DELTA>.

## Reproducibility Principle

Paper-facing workflows should read directly from public source files in `data/` whenever possible and avoid hidden `.RData` objects or manually generated intermediate TSV files. Serialized package data are reserved for documented external reference data that are part of the analysis provenance.


## Manuscript Trajectory Reproduction Notes

The package is designed to reproduce the manuscript trajectory analyses before generalizing them. Several details are intentionally manuscript-specific because they affect numerical equivalence with the original Figure 2 trajectory objects.

For the protein Figure 2C-style tumor versus normal/reference curves:

- tumor protein values are harmonized to public clinical sample IDs and collapsed to gene-level features by row mean when multiple protein rows map to the same gene;
- tumor matrices are row-scaled, then re-centered and re-scaled using tumor samples in the manuscript centering age range, `0 <= age <= 50`;
- trajectory fitting follows the original `proteo_tadj50.R` behavior and uses base `stats::loess()` defaults rather than `loess.control(surface = "direct")`;
- the original protein trajectory run fit tumor curves using samples with `age <= 80` (`age.cut2 = 80`) and excluded the public age-88 sample `C3N-02255` from the trajectory fit;
- Figure 2C curves are evaluated at all tumor sample ages `<= 50`, not only at ages from the sex being modeled. Male and female models are sex-stratified, but both are predicted on the same age support;
- normal/reference protein values are adjusted as in the manuscript code with `score ~ pH + PMI + Ethnicity`; missing normal pH is imputed once using the full normal-reference pH mean before sex stratification; ethnicity label `H` is combined with `C`, and `C` is used as the reference level;
- feature/sex/tissue-specific adaptive loess spans are part of the manuscript trajectory object and should be supplied explicitly when reproducing a published panel;
- confidence intervals are computed from the loess standard errors as `fit +/- qnorm(0.975) * se`.

These choices were checked against the original `tn_df.RData` and `protein_tadj50_list.RData` objects for `CNTN1`, `MAPT`, and `L1CAM`. With the manuscript settings above, package-generated Figure 2C `fit`, `low`, and `hi` values match the legacy `tn.df$value`, `tn.df$low`, and `tn.df$hi` to floating-point tolerance.
