# foundin_qtl
FOUNDIN-PD multi-omics QTL analysis

## This project contains the per-omic quantitative trait loci (QTL) analyses
- rnab-qtl: bulk RNAseq eQTL analysis
    - includes interaction eQTL analysis for bulk RNA with [SCADEN](https://github.com/KevinMenden/scaden) estimated cell-type fractions
    - includes quantitative trait score (QTS) analysis based on PD risk
    - includes colocalization with common PD risk loci
- atac-qtl: bulk ATACseq QTL (caQTL) analysis
    - includes quantitative trait score (QTS) analysis based on PD risk
    - includes colocalization with common PD risk loci
- scrn-qtl: single-cell RNAseq eQTL analysis
    - includes quantitative trait score (QTS) analysis based on PD risk
    - includes colocalization with common PD risk loci
- rnab3a-qtl: alternative polyadenylation (APA) quantative trait loci (3'aQTL) analysis from bulk RNAseq
    - includes quantitative trait score (QTS) analysis based on PD risk
    - includes colocalization with common PD risk loci
    
## refactor to remove or reduce duplicate code between different omic QTL analyses
### format modality sample info 
this may have already be done as part of another analysis project
1. combine known subject, sample, and modality covariates into single info table; quantifications/format_quants_covariates.ipynb
2. add data split column to the modality sample info table where the same test set subjects are used regardless of modality and day, test set based on batch 1 that was specifically balanced for the study, other samples are markered as training or exclude for known exclusions; analyses/format_data_splits_modeling.ipynb
Note: PDUI and RNAS are processed from the RNAB data so their info files are the same at start

### prepare the quantified modality
1. Reformat modalities input files as neccessary so that each modality can use same/similar prep and analysis code.
    - For PDUI reformat PDUI input from DaPar2 into quantified matrix and feature annotations, also correct the RNAB naming as needed; quantifications/dapar2_to_quant.ipynb
2. Scale and covariate adjust each modality across days without adjusting for differentiation state (day effect) and split into individual day data before and after scaling and adjustment; quantifications/split_quants_by_day.ipynb. Each modality can be run using quantifications/pm_run_split_quants_by_day.ipynb
3. Within each day scale and adjust each modalify within a specific day; quantifications/prep_quants_by_day.ipynb and quantifications/pm_run_prep_quants_by_day.ipynb

### <i>cis</i>-QTL analysis of prepared modality
1. Format the genotypes, current version of code was based on tensorQTL that required Plink bfile format, the newer version of tensorQTL may work with vcf but these analysis notebooks still expect bfiles. This only has to be performed once for all modalities and does not require removing unmatched modality samples, that is checked and performed in the <i>cis</i>-QTL notebook. The format genotypes notebook is genotypes/frmt_tensorqtl_genos.ipynb.
2. Run the <i>cis</i>-QTL analysis for a modality and day. This notebook finalized prep of inputs to ensure matched samples between genotype, modality, and specified covariates. The notebook runs tensorQTL <i>cis</i> map_nominal, map_cis with empirical, and map_independent; analyses/cis_tensorqtl.ipynb. This notebook can be run for all iterations of day and modalities with Papermill; analyses/pm_run_cis_tensorqtl.ipynb
3. Run the <i>cis</i>-interaction-QTL for the bulk modalities using the DA neuron fraction, NEEDS to be done

### Post processing of <i>cis</i>-QTL results
1. Compare <i>cis</i>-QTL between differentiation days for each modality; 
2. Scan <i>cis</i>-QTL results for intersection with Parksinson's disease risk for each differentiation day and modality
3. Colocalization analysis between Parksinson's disease risk and each differentiation day and modality <i>cis</i>-QTL results
4. Visual colocalization results

### Meta SCRN <i>cis</i>-eQTL analysis for FOUNDIN-PD (day 65) and Jerber et. al. HIPSCI (day 52) DA neuron single-cell <i>cis</i>-eQTL results
1. Prepare summary-stats for inclusion in meta-analysis
2. meta-analysis
3. Colocalization analysis between Parksinson's disease risk and final differentiation day(s) DA neuron <i>cis</i>-QTL results



