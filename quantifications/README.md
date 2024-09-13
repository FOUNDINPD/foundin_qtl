# prepare the quantified modality
1. combine known subject, sample, and modality covariates into single info table; format_quants_covariates.ipynb
2. Reformat modalities input files as neccessary so that each modality can use same/similar prep and analysis code.
    - ATAC use peaks_to_quant.ipynb
    - RNAB use edgeR_cpm_to_quant.ipynb
    - PDUI reformat PDUI input from DaPar2 into quantified matrix and feature annotations, also correct the RNAB naming as needed; dapar2_to_quant.ipynb
    - RNAS use exceRpt_to_quant.ipynb
    - CIRC use ciri2_to_quant.ipynb
    - METH use meffil_to_quant.ipynb
    - SCRN data already formated from different project
3. Split bulk modalities by day and single-cell modality by cell-type.
    - Bulk data, scale and covariate adjust each modality across days without adjusting for differentiation state (day effect) and split into individual day data before and after scaling and adjustment; split_quants_by_day.ipynb. Each modality can be run using pm_run_split_quants_by_day.ipynb
    - Single-cell use scrn_quants_by_celltype.ipynb
4. Within each day and modality scale and adjust by covariates generated from variance within the modality and day; prep_quants_by_day.ipynb and pm_run_prep_quants_by_day.ipynb