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

