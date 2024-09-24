# <i>cis</i>-QTL analysis of prepared modality
1. Format the genotypes, current version of code was based on tensorQTL that required Plink bfile format, the newer version of tensorQTL may work with vcf but these analysis notebooks still expect bfiles. This only has to be performed once for all modalities and does not require removing unmatched modality samples, that is checked and performed in the <i>cis</i>-QTL notebook. The format genotypes notebook is genotypes/frmt_tensorqtl_genos.ipynb.
2. Run the <i>cis</i>-QTL analysis for a modality and day. This notebook finalized prep of inputs to ensure matched samples between genotype, modality, and specified covariates. The notebook runs tensorQTL <i>cis</i> map_nominal and map_cis with empirical; analyses/cis_tensorqtl.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_cis_tensorqtl.ipynb
3. Run the <i>cis</i>-interaction-QTL for the bulk modalities using the DA neuron fraction, NEEDS to be done
4. Compare <i>cis</i>-QTL between differentiation days for each modality; analyses/compare_day_qtl_results.ipynb. This notebook can be run for all modalities with multi-day data using Papermill; analyses/pm_run_compare_day_qtl.ipynb

# Meta SCRN <i>cis</i>-eQTL analysis for FOUNDIN-PD (day 65) and Jerber et. al. HIPSCI (day 52) DA neuron single-cell <i>cis</i>-eQTL results
- Prepare summary-stats for inclusion in meta-analysis and run; analyses/meta_qtl.ipynb



# Colocalization of <i>cis</i>-QTL with common Parkinson's disease risk.
Using risk summary statistics from Nalls et al for PD (including 23andMe)
1. Compute LD statictics between risk index variant and other variants on the sample chromosome using AMP-PD WGS TOPMed freeze9 genotype calls as a reference panel. ld_risk_index_variant.ipynb
2. Identify ATAC peaks present in da65 bulk ATAC and the SCAT DA neuron cell-type that contain possible risk variants associated with PD. peak_risk_intersect.ipynb
1. Scan <i>cis</i>-QTL results for intersection with Parksinson's disease risk for each differentiation day and modality results; analyses/qtl_scan_risk.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_qtl_scan_risk.ipynb
2. Colocalization analysis between Parksinson's disease risk and each differentiation day and modality <i>cis</i>-QTL results including the meta-DAn eQTL results; analyses/colocalization.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_colocalization.ipynb
3. 
4. Visual colocalization results; figures/colocalization_results_heatmap.ipynb
    - includes Bryois et al single-cell eQTL from cortex cell types; analyses/coloc_bryois_brain_eqtl_results.ipynb
    - includes the meta-eQTL results of DA neurons and (p)bulk RNA (see below)

2. Colocalization analysis between Parksinson's disease risk and final differentiation day(s) DA neuron <i>cis</i>-QTL results; analyses/coloc_meta_eqtl_results.ipynb
    - also meta-analysis of FOUNDIN-PD da65 RNAB and HipSci D52 pseudobulk eQTL results
  
# <i>cis</i> proximal correlations between ATAC peaks and gene expression; Gene ~ ATAC, Gene ~ DNA methylation
1. Run the linear regressions between genes and their <i>cis</i> proximal ATAC peaks and genes and there <i>cis</i> proximal DNA methylation CpG sites. Here also using tensorQTL where ATAC peak features or CpG sites are the exogenious variable instead of variant genotype. cis_correlation.ipynb and cis_correlation_runner.ipynb
2. Compare results between cell-types. compare_celltype_ciscorr_results.ipynb
3. Identify <i>cis</i> ATAC peaks that correlated with each other by cell-type using Cicero. Monocle3_Cicero.ipynb
    
# QTS analysis
- Compute the genetic risk score for each of the FOUNDIN-PD cohort subjects based on the PD GWAS risk index variants or a subset of those variants that possible coloclize with DAn-meta eQTL results; genotypes/calculate_grs.ipynb
- Perform the QTS analysis, linear regression between GRS and modalities quantified features; analyses/day_qts.ipynb. Run this for each differentiation day and modality with Papermill; analyses/pm_run_day_qts.ipynb
- compare between days and modalities, none; analyses/compare_day_qts_results.ipynb
- visual QTS for specific feature; figures/plot_qts_pair.ipynb