# foundin_qtl
FOUNDIN-PD multi-omics QTL analysis

## This project contains the per-omic quantitative trait loci (QTL) analyses
- ATAC: caQTL for da0, da25, and da65
- METH: mQTL for da0 and da65
- RNAB: eQTL, 3'aQTL (PDUI) and circRNA eQTL (CIRC) for da0, da25, and da65 
- RNAS: eQTL (for miRNA, piRNA, and tRNA) for da0, da25, and da65 
- SCRN: eQTL for DA, ElC, eNP, iDA, lNP, NlC and PFPP cell-types and 3'aQTL (PDUI) for DA and iDA cell-types at da65

## This project also contains the per-omic quantitative trait score (QTS) analyses for each of the modality above as well
    
## Note: This code was refactored and consolidated to remove or reduce duplicate code between different omic QTL analyses after completion and the FOUNDIN-PD resource paper and the current analyses for the FOUNDIN-PD QTL project. The resource paper code relative to this repo is tagged as v0.1.0-alpha.

### prepare the quantified modality; [quantifications/](https://raw.githubusercontent.com/FOUNDINPD/foundin_qtl/quantifications/main/README.md)
1. combine known subject, sample, and modality covariates into single info table; format_quants_covariates.ipynb
2. Reformat modalities input files as neccessary so that each modality can use same/similar prep and analysis code.
    - ATAC use peaks_to_quant.ipynb
    - RNAB use edgeR_cpm_to_quant.ipynb
    - PDUI reformat PDUI input from DaPar2 into quantified matrix and feature annotations, also correct the RNAB naming as needed; dapar2_to_quant.ipynb
    - RNAS use exceRpt_to_quant.ipynb
    - CIRC use ciri2_to_quant.ipynb
    - METH data already formated from different project
    - SCRN data already formated from different project
3. Split bulk modalities by day and single-cell modality by cell-type.
    - Bulk data, scale and covariate adjust each modality across days without adjusting for differentiation state (day effect) and split into individual day data before and after scaling and adjustment; split_quants_by_day.ipynb. Each modality can be run using pm_run_split_quants_by_day.ipynb
    - Single-cell use scrn_quants_by_celltype.ipynb
4. Within each day and modality scale and adjust by covariates generated from variance within the modality and day; prep_quants_by_day.ipynb and pm_run_prep_quants_by_day.ipynb

### Analyses; [quantifications/](https://raw.githubusercontent.com/FOUNDINPD/foundin_qtl/analyses/main/README.md)
#### <i>cis</i>-QTL analysis of prepared modality
1. Format the genotypes, current version of code was based on tensorQTL that required Plink bfile format, the newer version of tensorQTL may work with vcf but these analysis notebooks still expect bfiles. This only has to be performed once for all modalities and does not require removing unmatched modality samples, that is checked and performed in the <i>cis</i>-QTL notebook. The format genotypes notebook is genotypes/frmt_tensorqtl_genos.ipynb.
2. Run the <i>cis</i>-QTL analysis for a modality and day. This notebook finalized prep of inputs to ensure matched samples between genotype, modality, and specified covariates. The notebook runs tensorQTL <i>cis</i> map_nominal and map_cis with empirical; analyses/cis_tensorqtl.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_cis_tensorqtl.ipynb
3. Run the <i>cis</i>-interaction-QTL for the bulk modalities using the DA neuron fraction, NEEDS to be done
4. Compare <i>cis</i>-QTL between differentiation days for each modality; analyses/compare_day_qtl_results.ipynb. This notebook can be run for all modalities with multi-day data using Papermill; analyses/pm_run_compare_day_qt.ipynb
5. For the single-cell eQTL compare eQTL between cell-types; analyses/compare_cell_qtl_results.ipynb

#### Meta SCRN <i>cis</i>-eQTL analysis for FOUNDIN-PD (day 65) and Jerber et. al. HIPSCI (day 52) DA neuron single-cell <i>cis</i>-eQTL results
- Prepare summary-stats for inclusion in meta-analysis and run meta-DAn eQTL; analyses/meta_qtl.ipynb


#### Colocalization of <i>cis</i>-QTL with common Parkinson's disease risk.
1. Format the meta-GWAS PD summary stats for colocaliation analysis; analyses/format_pd_sumstats.ipynb
1. Scan <i>cis</i>-QTL results for intersection with Parksinson's disease risk for each differentiation day and modality results; analyses/qtl_scan_risk.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_qtl_scan_risk.ipynb
2. Colocalization analysis between Parksinson's disease risk and each differentiation day and modality <i>cis</i>-QTL results where modalities include single-cell eQTL results as well; analyses/colocalization.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_colocalization.ipynb
3. Colocalization analysis between Parksinson's disease risk and meta-DAn eQTL; analyses/coloc_meta_eqtl_results.ipynb
4. Colocalization analysis between Parksinson's disease risk and public eQTL:
    - Bryois et al single-cell eQTL; analyses/coloc_bryois_brain_eqtl_results.ipynb
    - Jerber et al HiPSCI differentiated DAn eQTL; analyses/coloc_hipsci_dan_eqtl_results.ipynb
    - de Klein et al Metabrain bulk brain region eQTL; analyses/coloc_metabrain_eqtl_results.ipynb
  
#### <i>cis</i> proximal correlations between bulk modalities, such as Gene ~ ATAC, Gene ~ DNA methylation
1. Run the linear regressions between modality features and their <i>cis</i> proximal features from the other bulk modalities. Here also using tensorQTL where other modalities are the exogenious variable instead of variant genotype. cis_correlation.ipynb and pm_run_cis_correlation.ipynb
2. Compare results between differentiation days. compare_day_ciscorr_results.ipynb and pm_run_compare_days_ciscorr.ipynb
4. Compute LD statictics between risk index variant and other variants on the sample chromosome using AMP-PD WGS TOPMed freeze9 genotype calls as a reference panel. ld_risk_index_variant.ipynb
5. Identify ATAC peaks present in da65 bulk ATAC and the SCAT DA neuron cell-type that contain possible risk variants associated with PD. peak_risk_intersect.ipynb
6. Find the modality features that have a corellated <i>cis</i> proximal ATAC peak(s) where the peaks contain possible PD risk variants; ciscorr_scan_risk.ipynb and pm_run_ciscorr_scan_risk.ipynb


#### QTS analysis
- Compute the genetic risk score for each of the FOUNDIN-PD cohort subjects based on the PD GWAS risk index variants or a subset of those variants that possible coloclize with DAn-meta eQTL results; genotypes/calculate_grs.ipynb
- Perform the QTS analysis, linear regression between GRS and modalities quantified features; analyses/day_qts.ipynb. Run this for each differentiation day and modality with Papermill; analyses/pm_run_day_qts.ipynb
- compare between days and modalities, none; analyses/compare_day_qts_results.ipynb and analyses/pm_run_compare_day_qt.ipynb

#### Other analyses
- Features associated with differentiation day or estimated cell-type fraction. Run mixed effects model for each modality with subject as the random factor using samples from all timepoints the modality features and endogenous term and day or estimated DAn fraction as the exogenous term; mixed_effects_analysis_feature.ipynb and pm_run_mixed_effects.ipynb
- Similarly for known monogenic genes run mixed effects model where measures of cell-type differentiation are the endogenous term and the monogenic genes are exogenous terms; cell_fracs_corr_cis_monogenics.ipynb and pm_run_cell_fracs_cis_monogenics.ipynb

### Visualizations
