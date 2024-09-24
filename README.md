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

### <i>cis</i>-QTL analysis of prepared modality
1. Format the genotypes, current version of code was based on tensorQTL that required Plink bfile format, the newer version of tensorQTL may work with vcf but these analysis notebooks still expect bfiles. This only has to be performed once for all modalities and does not require removing unmatched modality samples, that is checked and performed in the <i>cis</i>-QTL notebook. The format genotypes notebook is genotypes/frmt_tensorqtl_genos.ipynb.
2. Run the <i>cis</i>-QTL analysis for a modality and day. This notebook finalized prep of inputs to ensure matched samples between genotype, modality, and specified covariates. The notebook runs tensorQTL <i>cis</i> map_nominal and map_cis with empirical; analyses/cis_tensorqtl.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_cis_tensorqtl.ipynb
3. Run the <i>cis</i>-interaction-QTL for the bulk modalities using the DA neuron fraction, NEEDS to be done

### Post processing of <i>cis</i>-QTL results
1. Compare <i>cis</i>-QTL between differentiation days for each modality; analyses/compare_day_qtl_results.ipynb. This notebook can be run for all modalities with multi-day data using Papermill; analyses/pm_run_compare_day_qtl.ipynb
2. Scan <i>cis</i>-QTL results for intersection with Parksinson's disease risk for each differentiation day and modality; analyses/qtl_scan_risk.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_qtl_scan_risk.ipynb
3. Colocalization analysis between Parksinson's disease risk and each differentiation day and modality <i>cis</i>-QTL results; analyses/colocalization.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_colocalization.ipynb
4. Visual colocalization results; figures/colocalization_results_heatmap.ipynb
    - includes Bryois et al single-cell eQTL from cortex cell types; analyses/coloc_bryois_brain_eqtl_results.ipynb
    - includes the meta-eQTL results of DA neurons and (p)bulk RNA (see below)

### Meta SCRN <i>cis</i>-eQTL analysis for FOUNDIN-PD (day 65) and Jerber et. al. HIPSCI (day 52) DA neuron single-cell <i>cis</i>-eQTL results
1. Prepare summary-stats for inclusion in meta-analysis and run; analyses/meta_qtl.ipynb
2. Colocalization analysis between Parksinson's disease risk and final differentiation day(s) DA neuron <i>cis</i>-QTL results; analyses/coloc_meta_eqtl_results.ipynb
    - also meta-analysis of FOUNDIN-PD da65 RNAB and HipSci D52 pseudobulk eQTL results
    
### QTS analysis
- Compute the genetic risk score for each of the FOUNDIN-PD cohort subjects based on the PD GWAS risk index variants or a subset of those variants that possible coloclize with DAn-meta eQTL results; genotypes/calculate_grs.ipynb
- Perform the QTS analysis, linear regression between GRS and modalities quantified features; analyses/day_qts.ipynb. Run this for each differentiation day and modality with Papermill; analyses/pm_run_day_qts.ipynb
- compare between days and modalities, none; analyses/compare_day_qts_results.ipynb
- visual QTS for specific feature; figures/plot_qts_pair.ipynb



