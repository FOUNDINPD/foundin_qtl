These are the steps for running a simple QTL analysis on the ATAC data where each 'day' is run separately; ie simple style.

These analyses uses a JupyterLab instance(s) to run the notebook set on a local instance or Google Cloud Platform (GCP) AI Platform Notebooks instance with decent resources. Such as n1-standard-32 (32 vCPUs, 120 GB RAM) with NVIDIA Tesla T4 GPU, usually with 500GB or 1TB persistent disk. 

Using tensorQTL module to run QTL, based on pytorch.

The GCP cloud AI notebooks interface is here [AI Platform Notebooks](https://console.cloud.google.com/ai-platform/notebooks/list/instances?project=foundin-pd). Select your appropriate project. The GCP documentation for here [AI Platform Notebooks documentation](https://cloud.google.com/ai-platform/notebooks/docs).

If using GCP AI notebooks once the instance has been created you should see the OPEN JUPYTERLAB link to access the JupyterLab instance.

Within JupyterLab open a Terminal
1. pull the code from github [FOUNDIN-PD caQTL](foundin-pd_caqtl)
    - recommend keeping code dir and data/analyses directories separate
    - note below the 'notebooks' directory is the github cloned repo directory

2. open the notebooks/setup_and_data_pull.ipynb and run, this will create the analysis data directories and pull input data such as genotypes, quantified features, feature annotations, and some sample info. Have to set the appropriate GCP bucket paths; pointing to source data and specify the root analysis directory that will be created

3. run notebooks/format_quants_covariates.ipynb, which will pull the various known subject and sample covariates, clean them up as needed and save into a single covariate file, appropriate for this assay's analysis, currently using:
     - AMP-PD subject info: amppd_demographicsPlus_2019_v1release_1015.csv
     - FOUNDIN-PD cell line info: cell_metadata.csv
     - Genetics PCs generated from PPMI WGS genotypes: foundin.freeze9.pca.eigenvec
     - estimated cell fractions from RNAB using SCADEN deconvolution: rnab_cell_fracs_scaden.csv
     - atac metrics: foundin_atac_metrics.csv
     - FOUNDIN-PD subject info: Expanded_overview_of_included_PPMI_samples_overview.csv
     - FOUNDIN-PD PD risk GRS: Expanded_overview_of_included_PPMI_samples_GRS.csv

4. run the notebooks/split_quants_by_day.ipynb, which will split the quantified data by days present in the data, ie for ATAC should be da0, da25, and da65. This notebooks also does covariate generation from variance within the quantified features and compare with known sample covariates for data across the differentiation timepoints.

5. run notebooks/frmt_tensorqt_genos.ipynb, tensorqtl uses panda_plink to load the genotypes which basically reads plink bfiles as dosages, since we typically just use per chrom vcfs or plink pfiles need to convert to plink pfiles, and since cohort is pretty small and tensorQTL pretty fast go ahead and merge the per chrom pfiles into autosome bfile
    
- note since most analysis and prep are repeated over the three separate cell differentiaion timepoints individually I used [Papermill](https://papermill.readthedocs.io/en/latest/) to automate the notebook runs instead of manually creating three separate experiment notebooks. So you will have the 'template' notebooks pulled from github repo but then for analyses repeated by day those will have generated notebooks, which will be found in notebooks/pm_gend_nbs

6. Prep the quantified feature data by running notebooks/quants_prep_runner.ipynb. This notebook uses papermill to run notebooks/day_quants_prep.ipynb where it just changes the day variable (parameter), it will then create an output notebook for that run (experiment). This notebook only work with samples for the specified differentiation day. It subsets the quantified features based on detection rate (here 75%), quantile transform and scales the data, and adjusts the quantified features by covariates generated from the quantified features using UMAP projections from the quantified features. The UMAP projection is based on features after excluding features from the bottom quartile of variance and those that are 'detected'. The notebook formats the  quantified and scaled features matrix as a Hierarchical Data Format version 5 (HDF5) files and then also phenotype bed file for use by tensorQTL. 

    - Note that in the above prep step since the generated and adjusted quantified data is based on covariates generated from the data, this is one reason why the individual differention days are prepped separately. As would we expected from differntiating cells the timepoint and cell-type content would be the primary sources of variation if the data was done together. The goal is to remove confounding known and unknown sources of variation for the data not the interesting differentiated cell-specific variation.

    - Review the quantified data using the Papermill generated notebooks. For instance the notebook generated and run for the d25 quantified data would be found notebooks/pm_gend_nbs/foundin_d25_quantified_prep.ipynb. Here everthing should be fine but prior to running the QTL analysis would want to verify how things looked as far as how many features detected, scaling transformation doesn't look off, and how the adjustment for generated covariates changed the data (also includes some information on what known covariates are correlated with the generated ones). Check to see if the generated covaritates are going to remove signal related to variables that may be used in other analyses for future reference; ie disease status. Here for simple QTL where detecting the effect of genotype on quantified feature it is OK if covariate removes that signal. So, ideally looking at closer to pure genotype effect.

7. Run the simple QTL analysis, again using Papermill to run the days separately from same template notebook, by running the notebooks/day_cis_qtl_runner.ipynb notebook. The primary outputs here are the generated notebooks, such as notebooks/pm_gend_nbs/day_cis_qtl_tensorqtl.ipynb, all analysis results per chromosome as parquet files, the top result per feature, and all significant independent results based on FDR corrections implemented in tensorQTL. 

8. compare QTL results by differentiation day
    notebooks/compare_day_indep_results.ipynb
    
9. scan QTL results for PD meta5 risk variants
    notebooks/eqtl_scan_risk.ipynb
    
10. for any PD risk and QTL colocalization plot those specific features
    notebooks/feature_specific_qtl_results_runner.ipynb
    
11. plot specific QTL variants
    notebooks/plot_qtl_genotype.ipynb
    
12. run QTS analysis using PD GRS based on meta5 results
    foundin/notebooks/calculate_grs.ipynb
    foundin/notebooks/day_qts_runner.ipynb
    foundin/notebooks/compare_day_qts_results.ipynb

13. run longitudinal feature analysis, this runs linear mixed effects model to find features over differentiation days where subject is the random effect; ie feature ~ day + (1|subject) also runs mixed effects model plus couple covariates for cell fraction and UMAPs covariates
    notebooks/longitudinal_analysis_feature.ipynb
    
14. To button everything up run the notebooks/finish_and_push_back.ipynb notebook to push results, info files, and generated notebooks back to the Google cloud storage (GCS) buckets

15. If completely done you used a GCP instance should delete the compute instance from the GCP AI Platform Notebooks interface. Note if you feel like you need to do some more exploratory analysis or poking around you can just stop the instance instead of deleting and then come back later a start the instance again and jump back in. The instance disks are persistent, obviously memmory state isn't presevered, but data and tooling on the disk are. However while the cost is much reduced by not having the compute resources running you are still paying small fee in storage for the persistent disk and machine image so don't leave in place if not needed any longer. By the way sometimes helpful, you can stop the instance and change the machine type resources and restart so if need more or less cpu, memory, even gpu can do.


