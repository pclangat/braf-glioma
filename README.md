# Integrated molecular and clinical analysis of BRAF-mutant glioma in adults 

This repository holds the analysis and data visualization codes used in:

Schreck KC, Langat P, Bhave VM, Li T, Woodward E, Pratilas CA, Eberhart CG, Bi WB (2022). Integrated molecular and clinical analysis of BRAF-mutant glioma in adults . npj Precision Oncology, in press.


## Summary of analyses by Figure
**Figure 1**
- Age cutpoint analysis: `braf_age_cutpoints.R`
- Clinical data: `braf_mastersheet_compiler.R` --> `braf_clinical_composite.R`



**Figure 2**
- Genomic characterization, TMB calculation, and oncoprint generation: `braf_oncoprint.R` and `MafFunctions_BRAF.R`



**Figure 3**
- RNAseq analysis: `braf_RNAseq.Rmd` --> `braf_logisticPCA_kmeans_clustering.Rmd` --> `braf_heatmap.Rmd`
- Correlation testing: 'braf_chisq_correlations.R'
- Covariate testing: `braf_covariates_analysis.R`
- Fusions circos: `braf_circos.R`



**Figure 4**
- Survival analysis: `braf_survival.R`



**Figure 5**
- Swimmers plot: `braf_swimmers.R`
