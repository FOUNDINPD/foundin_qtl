"""
    Utility functions for use in notebooks related to prepping the modality quantifications
"""

# imports
from pandas import DataFrame , get_dummies
from numpy import cumsum
from sklearn.preprocessing import QuantileTransformer, MinMaxScaler
from seaborn import distplot , scatterplot, heatmap
from random import sample
from umap import UMAP
import ppscore as pps
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from scipy.stats import zscore
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

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
    functions to generate and visualize known and unknow covariates using PCA, UMAP and PPScore
"""
# function for plotting umap of traits with covar high lights
def plot_umap_clusters(umap_df, hue_cov=None, style_cov=None, size_cov=None):
    with rc_context({'figure.figsize': (12, 12), 'figure.dpi': dpi_value}):
        plt.style.use('seaborn-bright')
        sns_plot = scatterplot(x='x_umap',y='y_umap',
                               hue=hue_cov, style=style_cov, size=size_cov,
                               data=umap_df)
        plt.xlabel('x-umap')
        plt.ylabel('y-umap')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0,prop={'size': 10})
        plt.show()

def plot_pca_pair(this_df: DataFrame, first: str, second: str, 
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
        
# small function to generate umap from pandas dataframe, for all features (columns)
# and return back as dataframe with source index intact
def generate_umap_covs_df(this_df, other_covs_df=None, 
                             rnd_digits=3, merge_input=False):
    #run UMAP on the data frame features
    umap_results = UMAP(random_state=42).fit_transform(this_df)
    umap_df = DataFrame(umap_results,columns=['x_umap','y_umap'],
                                       index=this_df.index).round(rnd_digits)
    if merge_input:
        umap_df = umap_df.merge(this_df,left_index=True,right_index=True)
    if other_covs_df is not None:
        umap_df = umap_df.merge(other_covs_df, how='left', 
                                left_index=True, right_index=True)
    print(f'The dimensions of the umap df and the traits are {umap_df.shape}')
    return umap_df 

def generate_pca(this_df: DataFrame, plot_variance: bool=False, 
                 amount_var: float=0.99, rnd_digits: int=3) -> DataFrame:
    pca = PCA(n_components=amount_var)
    components = pca.fit_transform(this_df)
    pca_df = DataFrame(data=components, index=this_df.index).round(rnd_digits)
    pca_df = pca_df.add_prefix('PC_')
    print(pca.explained_variance_ratio_)
    print(len(pca.explained_variance_ratio_))    
    if plot_variance:
        with rc_context({'figure.figsize': (8, 8), 'figure.dpi': dpi_value}):
            plt.style.use('seaborn-bright')
            plt.plot(cumsum(pca.explained_variance_ratio_))
            plt.xlabel('Number of components')
            plt.ylabel('Explained variance')
            plt.show()
    return pca_df

# function to iterate over target features and use PPScore to find covarites of interest
def pps_predict_targets(this_df, target_list, min_ppscore: float=0.05):
    covs_to_check = []
#     covs_list = ['x_umap', 'y_umap']
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
    matrix_df = pps.matrix(this_df[(set(covs_to_check) | set(cov_targets))])
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
def plot_correlation_heatmap(this_df, covs_list : list=None, min_pearson: float=0.22):
    cor = this_df.corr(method='pearson')
    cor.dropna(how='all', inplace=True)
    modified_title = ''
    if covs_list is not None:
        
        limited_cor = cor[covs_list]
        cor = limited_cor.loc[(limited_cor['x_umap'].abs() > min_pearson) | 
                              (limited_cor['y_umap'].abs() > min_pearson)]
        modified_title = 'limited'
    print(cor.shape)
    fig_width = cor.shape[1] if cor.shape[1] > 12 else 12
    fig_height = cor.shape[0] if cor.shape[1] > 12 else 12
    with rc_context({'figure.figsize': (fig_width, fig_height), 'figure.dpi': dpi_value}):
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
    dums_df = get_dummies(cats_df)
    print(f'one-hot encoded categoricals shape {dums_df.shape}')

    temp_df = temp_df.merge(dums_df, how='inner', left_index=True, right_index=True)
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
def scale_dataframe(this_df : DataFrame):
    scaledX = MinMaxScaler().fit_transform(QuantileTransformer(output_distribution='normal')
                                           .fit_transform(this_df))
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
    quants_wd_var_df = this_df[keep_ids]
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