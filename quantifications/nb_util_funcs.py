"""
    Utility functions for use in notebooks related to prepping the modality quantifications
"""

# imports
from pandas import DataFrame , get_dummies
from numpy import cumsum, arange, where
from sklearn.preprocessing import QuantileTransformer, MinMaxScaler
from seaborn import distplot, scatterplot, heatmap
from random import sample
from umap import UMAP
import ppscore as pps
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA, FastICA, NMF
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import zscore
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from kneed import KneeLocator
import pymde
import torch

# variables
dpi_value = 100

# Notebook functions
"""
    functions for detection rates calculations and plotting
"""

def calculate_detection_rates(this_df, modality, round_percision=1, 
                              min_quant_value=None):
    if min_quant_value is None:
        min_quant_value = this_df.round(round_percision).min().min()

    print(f'minimun {modality} value is {min_quant_value}')

    detected_df = this_df.mask(this_df.round(round_percision) <= min_quant_value, 0)

    # calculate the missing counts from the detected df mask
    trait_missing_rates = round(detected_df.isin({0}).sum(0)/detected_df.shape[0], 2)
    sample_missing_rates = round(detected_df.isin({0}).sum(1)/detected_df.shape[1], 2)

    print(f'{len(trait_missing_rates)} features with mean missing \
rate = {trait_missing_rates.mean()}')
    print(f'{len(sample_missing_rates)} samples with mean missing \
rate = {sample_missing_rates.mean()}')
    return trait_missing_rates, sample_missing_rates

def plot_missing_rates(feature_rates, sample_rates):
    with rc_context({'figure.figsize': (12, 12), 'figure.dpi': dpi_value}):
        plt.style.use('seaborn-bright')
        plt.subplot(2, 2, 1)
        distplot(feature_rates.values)
        plt.title('Features missingness rates')
        plt.subplot(2, 2, 2)
        distplot(sample_rates.values)
        plt.title('Samples missingness rates')
        plt.show()
    
def bad_callrate_features(features_missing_rates, max_missing_rate):
    bad_call_rates = features_missing_rates[features_missing_rates > max_missing_rate]
    print(f'features with bad call rates shape {bad_call_rates.shape}, \
fraction of features with bad rates {bad_call_rates.shape[0]/features_missing_rates.shape[0]}')
    return bad_call_rates

def subset_well_detected_features(this_df, bad_call_rates):
    detected_traits = list(set(this_df.columns)-set(bad_call_rates.index))
    this_wd_df = this_df[detected_traits]
    print(f'shape of well detected quants {this_wd_df.shape}')
    return this_wd_df

"""
    functions to generate and visualize known and unknow covariates using PCA, NMF, ICA and PPScore
"""
def iterate_model_component_counts(max_count: int, data_df: DataFrame, 
                                   model_type: str=['PCA', 'NMF', 'ICA']) -> (list, list):
    r2_rets = []
    rmse_rets = []
    for comp_num in arange(1, max_count+1):    
        _,_,r2,rmse = generate_selected_model(comp_num, data_df, model_type)
        r2_rets.append(r2)
        rmse_rets.append(rmse)
    return r2_rets, rmse_rets

def generate_selected_model(n_comps: int, data_df: DataFrame,
                            model_type: str=['PCA', 'NMF', 'ICA']) -> (object, DataFrame, float, float):
    if model_type == 'PCA':
        model = PCA(n_components=n_comps, random_state=42)
    if model_type == 'NMF':
        model = NMF(n_components=n_comps, init='random', random_state=42, max_iter=500)
    if model_type == 'ICA':
        model = FastICA(n_components=n_comps, random_state=42)        
    components= model.fit_transform(data_df)
    recon_input = model.inverse_transform(components)
    r2 = r2_score(y_true=data_df, y_pred=recon_input)
    rmse = mean_squared_error(data_df, recon_input, squared=False)
    print(f'{model_type} with {n_comps} components accuracy is {r2:.4f}, RMSE is {rmse:.4f}')  
    ret_df = DataFrame(data=components, index=data_df.index).round(4)
    ret_df = ret_df.add_prefix(f'{model_type}_')
    return model, ret_df, r2, rmse

def component_from_max_curve(scores, label: str=['R2', 'RMSE', 'EVR']) -> int:
    if label == 'R2' or label == 'EVR':
        data_curve = 'concave'
        data_direction = 'increasing'
    if label == 'RMSE':
        data_curve = 'convex'
        data_direction = 'decreasing'        
    knee = KneeLocator(arange(1, len(scores)+1), scores, 
                       S=1.0, curve=data_curve, direction=data_direction)
    if knee.knee is None:
        num_comp = len(scores)
    else:
        num_comp = int(knee.knee) 
    print(f'best curve at knee {num_comp}')
    exp_value = scores[num_comp-1]
    print(f'best number of components is {num_comp} at {label} of {exp_value}')
    knee.plot_knee()
    plt.show()
    knee.plot_knee_normalized()
    plt.show()
    return num_comp

def plot_pair(this_df: DataFrame, first: str, second: str, 
              hue_cov=None, style_cov=None, size_cov=None):
    with rc_context({'figure.figsize': (8, 8), 'figure.dpi': dpi_value}):
        plt.style.use('seaborn-bright')
        sns_plot = scatterplot(x=first,y=second, 
                               hue=hue_cov, style=style_cov, size=size_cov, 
                               data=this_df, palette='viridis')
        plt.xlabel(first)
        plt.ylabel(second)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, 
                   borderaxespad=0,prop={'size': 10})
        plt.show() 
        
# small function to generate 2D embed from pandas dataframe, for all features (columns)
def generate_2d_embed_df(this_df, covs_df=None, embed_type: str='MDE', rnd_digits=3, merge_input=False):
    ret_df = None
    if embed_type == 'MDE':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        print(device)
        pymde.seed(42)
        mde = pymde.preserve_neighbors(this_df.to_numpy(), device=device, verbose=True)
        embedding = mde.embed(verbose=True)
        embed_results = embedding.detach().cpu().numpy()
    elif embed_type == 'UMAP':
        embed_results = UMAP(random_state=42).fit_transform(this_df)

    ret_df = DataFrame(embed_results, columns=[f'{embed_type}_1',f'{embed_type}_2'],
                       index=this_df.index).round(rnd_digits)
    if merge_input:
        ret_df = ret_df.merge(this_df,left_index=True, right_index=True)
    if covs_df is not None:
        ret_df = ret_df.merge(covs_df, how='left', left_index=True, right_index=True)
    print(f'The dimensions of the embed df and the covariates are {ret_df.shape}')
    return ret_df

def show_2d_embed(data_df: DataFrame, cov_df: DataFrame, type: str='MDE',
                  hue: str=None, style: str=None, size: str=None, 
                  verbose: bool=False):
    embd_df = generate_2d_embed_df(data_df, cov_df, embed_type=type)
    print(f'embd_df shape is {embd_df.shape}')
    if verbose:
        display(embd_df.head())
    plot_pair(embd_df, f'{type}_1', f'{type}_2', 
              hue_cov=hue, style_cov=style, size_cov=size) 

# function to iterate over target features and use PPScore to find covarites of interest
def pps_predict_targets(this_df, target_list, min_ppscore: float=0.05):
    covs_to_check = []
    for this_cov in target_list:
        print(this_cov)
        predictors_df = pps.predictors(this_df, this_cov)
        # drop anything that has ppscore of zero
        predictors_df = predictors_df.loc[predictors_df['ppscore'] > min_ppscore]
        display(predictors_df)
        covs_to_check.extend(list(predictors_df['x'].values))

    print(f'found {len(covs_to_check)} covariates that may preditct target covariates')    
    return covs_to_check

# plot ppscore matrix 
def plot_ppscore_matrix(this_df, covs_to_check, cov_targets, min_ppscore: float=0.05):
    matrix_df = pps.matrix(this_df[list(set(covs_to_check) | set(cov_targets))])
    matrix_df = matrix_df.loc[matrix_df['ppscore'] > min_ppscore]
    print(matrix_df.shape)

    matrix_df['ppscore'] = matrix_df['ppscore'].round(2)
    plot_matrix_df = matrix_df[['x', 'y', 'ppscore']].pivot(columns='x', index='y', values='ppscore')
    print(plot_matrix_df.shape)
    # display(plot_matrix_df)
    with rc_context({'figure.dpi': dpi_value}):
        plt.style.use('seaborn-bright')
        plt.figure(figsize=(plot_matrix_df.shape[0],plot_matrix_df.shape[1])) 
        heatmap(plot_matrix_df, vmin=0, vmax=1, cmap='Blues', linewidths=0.05, 
                annot=True, annot_kws={'fontsize':12})
        plt.title('PPScore heatmap')
        plt.show()
    
# plot heatmap of Pearson correlation matrix for PPScore covariates
def plot_correlation_heatmap(this_df, min_pearson: float=0.22):
    cor = this_df.corr(method='pearson')
    cor.dropna(how='all', inplace=True)
    modified_title = ''
    print(cor.shape)
    fig_width = cor.shape[1] if cor.shape[1] > 12 else 12
    fig_height = cor.shape[0] if cor.shape[1] > 12 else 12
    with rc_context({'figure.figsize': (fig_width, fig_height), 'figure.dpi': 100}):
        plt.style.use('seaborn-bright')    
        heatmap(cor[(cor > min_pearson) | (cor < -min_pearson)], annot=True, 
                annot_kws={"fontsize":10}, linewidths=0.05, cmap='Blues')    
        plt.title(f'Pearson heatmap of PPScore covariates {modified_title}')
        plt.show()

# function to one-hot encode the categorical covariates and merge with continuous ones    
def dummy_covs_as_needed(this_df):
    temp_df = this_df.copy()
    cats_df = temp_df.select_dtypes(include=['object'])
    print(f'categoricals shape {cats_df.shape}')
    dums_df = get_dummies(cats_df).astype(int)
    print(f'one-hot encoded categoricals shape {dums_df.shape}')
    temp_df = temp_df.drop(columns=cats_df.columns)
    temp_df = temp_df.merge(dums_df, how='inner', left_index=True, 
                            right_index=True)
    print(f'new covs df shape {temp_df.shape}')
    return temp_df

"""
    additional visualization functions
"""
# small function to plot before and after of transform based on named feature,
# or if a feature isn't specified then one pull at random
def plot_trnsfrm_effect_example(before_df, after_df, feature_id=None,
                                bf_label='quantile transformed', 
                                af_label='quantile transformed and covariate adjusted'):
    # if no feature ID provided get randome one
    if feature_id is None:
        feature_id = sample(list(after_df.columns), 1)[0]
    with rc_context({'figure.figsize': (8, 8), 'figure.dpi': dpi_value}):
        plt.style.use('seaborn-bright')    
        distplot(before_df[feature_id])
        plt.title(f'{feature_id} {bf_label}')
        plt.show()
        distplot(after_df[feature_id])
        plt.title(f'{feature_id} {af_label}')
        plt.show()
        scatterplot(x=before_df[feature_id], y=after_df[feature_id])
        plt.title(f'{feature_id}')
        plt.xlabel(f'{bf_label}')
        plt.ylabel(f'{af_label}')
        
"""
    analysis functions
"""        
# small function to perform the quantile transform and minmax scale on a pandas dataframe
def scale_dataframe(this_df : DataFrame, with_qt: bool=True):
    if with_qt:
        scaledX = MinMaxScaler().fit_transform(QuantileTransformer(output_distribution='normal')
                                               .fit_transform(this_df))
    else:
        scaledX = MinMaxScaler().fit_transform(this_df)    
    scaled_df = DataFrame(data=scaledX, columns=this_df.columns, 
                                 index=this_df.index)  
    return scaled_df    

# exclude low variance features from covariate generation
def exclude_low_var_features(this_df: DataFrame, quartile_to_drop: str ='25%', 
                             known_feature_to_drop=None):
    quants_vars = this_df.var() 
    print(quants_vars.describe())
    # drop features within the lower quartile of variance
    min_variance = quants_vars.describe()[quartile_to_drop]
    # min_variance = quants_vars.describe()['50%']
    keep = quants_vars[quants_vars > min_variance]
    if known_feature_to_drop is not None:
        keep_ids = set(keep.index) - set(known_feature_to_drop)
    else:
        keep_ids = set(keep.index)
    quants_wd_var_df = this_df[list(keep_ids)]
    print(f'shape of the features to keep {keep.shape}')
    print(f'shape of input features df {this_df.shape}')
    print(f'shape of variance features df {quants_wd_var_df.shape}')
    return quants_wd_var_df

# function to fit linear model to covariates and calculate the standardized residuals
def covariate_residuals(traits_df, covars_df):
    lm = LinearRegression(n_jobs=16)
    residuals_df = traits_df.copy()
    covar_scores_by_trait = {}

    for trait in traits_df:
            model = lm.fit(covars_df, traits_df[trait])
            covar_scores_by_trait[trait] = model.score(covars_df,traits_df[trait])
            model_predicted = model.predict(covars_df)
            residuals_df[trait] = zscore(traits_df[trait] - model_predicted)
            
    # scale the residuals
    residualsX = MinMaxScaler().fit_transform(residuals_df)
    residuals_df = DataFrame(data=residualsX, columns=traits_df.columns, 
                                index=traits_df.index)
    
    # grab the covariates model scores
    covar_scores_by_trait_df = DataFrame.from_dict(covar_scores_by_trait,
                                                      columns=['score'],
                                                      orient='index').round(3)
    covar_scores_by_trait_df.index.name = 'featureID'
    return residuals_df, covar_scores_by_trait_df

"""
    input output functions
"""
# small function to save parquet file
def write_df_to_parquet(this_df, file_name):
    this_df.to_parquet(file_name)
    
# small function to save hdf file
def write_df_to_hdf(this_df, file_name, key='quants', mode='w'):
    this_df.to_hdf(file_name, key=key, mode=mode)

# format dataframes to tensorQTL pheno bed files
def format_tensorqtl_bed(data_df: DataFrame, features_df: DataFrame, out_file: str, 
                         modality: str, id_map: DataFrame):
    # get feature annots for present features
    feature_present_df = features_df.loc[features_df['feature'].isin(data_df.columns)]
    print(f'features present shape {feature_present_df.shape}')
    # tensorQTL pheno bed is rows = features and columns = samples
    # where first four columns are chr, start, end, phenotype_id, then sample1 ... sampleN

    # create dict for renaming columns (samples) from assayid to geno_id
    sample_col_dict = id_map.set_index('assayid').to_dict()['sampleid']
    
    # transpose the data df from sample x feature to feature x sample
    tdata_df = data_df.transpose()
    
    # modify annots
    feature_present_df = feature_present_df.rename(columns={'chrom': 'chr'})
    feature_present_df['end'] = feature_present_df['start'] + 1
    print(f'features present shape {feature_present_df.shape}')
    feature_present_df = feature_present_df.drop_duplicates(subset=['feature'], 
                                                            keep='first', 
                                                            ignore_index=True)
    feature_present_df = feature_present_df.set_index('feature', drop=False)
    feature_present_df = feature_present_df.reindex(tdata_df.index)
    
    # insert the feature annots
    tdata_df.insert( 0, column='chr', value=feature_present_df['chr'])
    tdata_df.insert( 1, column='start', value=feature_present_df['start'])
    tdata_df.insert( 2, column='end', value=feature_present_df['end'])
    tdata_df.insert( 3, column='phenotype_id', value=feature_present_df['feature'])
    
    # if there are any genes that were in quants but not feature annots
    # remove these with missing positions
    tdata_df = tdata_df.loc[~tdata_df['chr'].isna()]
    # make the positions ints instead of floats
    tdata_df['start'] = tdata_df['start'].astype('int64')
    tdata_df['end'] = tdata_df['end'].astype('int64')

    # now rename sample ids in columns
    tdata_df = tdata_df.rename(columns=sample_col_dict)
    
    # tensorQTL pheno reader expects sorted bed
    tdata_df = tdata_df.sort_values(['chr', 'start', 'end'])
    
    # save the file
    tdata_df.to_csv(out_file, index=False, sep='\t', compression='gzip')
    print(f'tdata_df shape is {tdata_df.shape}')