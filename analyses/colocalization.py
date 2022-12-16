"""
Colocalization related functions
These are python versions of methods and functions based on the coloc R library 
    that were ported to python by tensorQTL and OpenTargets package authors
    and then adapted for colocalization analyses here

coloc library code here
    https://github.com/chr1swallace/coloc
    https://github.com/chr1swallace/coloc/blob/main/R/claudia.R
coloc references
    https://doi.org/10.1371/journal.pgen.1004383
    https://doi.org/10.1371/journal.pgen.1008720
    https://doi.org/10.1101/2021.02.23.432421
    
OpenTargets finemapping code here
    https://github.com/opentargets/genetics-finemapping
    https://github.com/opentargets/genetics-finemapping/blob/master/finemapping/credible_set.py
OpenTargers reference
    https://doi.org/10.1038/s41588-021-00945-5
    
tensorQTL code here
    https://github.com/broadinstitute/tensorqtl
    https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/coloc.py  
tensorQTL reference
    https://doi.org/10.1186/s13059-019-1836-7
"""

import numpy as np
from scipy.stats import norm
from pandas import DataFrame, Series

# from or based on opentargets finemapping code
def calc_abf(pval: float, maf: float, n: int, prop_cases=None):
    """ Caluclate Approximate Bayes Factor (Wakefield, 2009, Genet Epidemiol.).
        Based on code from coloc: https://github.com/chr1swallace/coloc
    Args:
        pval (float): GWAS p-value
        maf (float): Minor allele freq
        n (int): Sample size
        prop_cases (float or None): number of cases, if left blank will assume
            quantitative trait
    Returns:
        natural log(ABF)
    """
    # Assert/set types
    pval = float(pval)
    maf = float(maf)
    n = int(n)
    prop_cases = float(prop_cases) if prop_cases else None

    # Estimate variance for quant trait
    if prop_cases is None:
        sd_prior = 0.15
        v = var_data(maf, n)
    # Estimate var for cc study
    else:
        sd_prior = 0.2
        v = var_data_cc(maf, n, prop_cases)

    # Calculate Z-score
    z = np.absolute(norm.ppf(pval / 2))

    # Calc shrinkage factor: ratio of the prior variance to the total variance
    r = sd_prior**2 / (sd_prior**2 + v)

    # Approximate BF - ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (np.log(1 - r) + (r * z**2))

    return lABF

def var_data(maf: float, n: int):
    """ Calc variance of MLE of beta for quantitative trait, assuming var(y)=0
    """
    var = 1 / (2 * n * maf * (1 - maf))
    return var

def var_data_cc(maf: float, n: int, prop_cases: float):
    """ Calc variance of MLE of beta for case-control
    """
    var = 1 / (2 * n * maf * (1 - maf) * prop_cases * (1 - prop_cases))
    return var

def freq_to_maf(freq: float):
    """ Convert allele frequency to minor allele freq
    """
    return min(freq, 1-freq)

def credible_sets(df: DataFrame):
    """ Identify the 95 and 99% credible sets
    Args:
        df (pandas.DataFrame) must include posterior probability (PP) column
    """    
    # Calc cumulative sum of the posterior probabilities
    df['pp_cumsum'] = df.PP.cumsum()

    # Find 99% and 95% credible sets - this is horrible
    set_idx = df.pp_cumsum.gt(0.95).tolist().index(True)
    df['is95_credset'] = [1] * (set_idx + 1) + [0] * (df.shape[0] - (set_idx + 1))
    df['is95_credset'] = df.is95_credset.map({1:True, 0:False})
    set_idx = df.pp_cumsum.gt(0.99).tolist().index(True)
    df['is99_credset'] = [1] * (set_idx + 1) + [0] * (df.shape[0] - (set_idx + 1))
    df['is99_credset'] = df.is99_credset.map({1:True, 0:False})

# from or based on tensorqtl code
def compute_pp(labf: Series) -> Series:
    """ Compute the posterior probabilities
    Args:
        labf (pandas.Series) series of ABF's to compute the PP's on
    """     
    sum_lABF = logsumexp(labf)
    ret_values = (labf - sum_lABF).apply(np.exp)
    return ret_values

def logsumexp(x):
    """ Calculates the log of the sum of the exponentiated logs taking out the
        max, i.e. insuring that the sum is not Inf
    Args:
        x (array_like)
    Returns:
        Sum of exponentiated logs
    """    
    mmax = np.max(x)
    return mmax + np.log(np.sum(np.exp(x-mmax)))

def logdiff(x, y):
    """ calculates the log of the difference of the exponentiated
    logs taking out the max, i.e. insuring that the difference is not negative
    Args:
        x (array_like)
        y (array_like)
    """
    mmax = max(x, y)
    return mmax + np.log(np.exp(x - mmax) - np.exp(y - mmax))

def combine_abf(l1, l2, p1: float=1e-4, p2: float=1e-4, p12: float=1e-5):
    """ Perform the colocalization analysis
    Args:
        l1 (array_like) trait1's ABFs
        l2 (array_like) trait2's ABFS
        p1 (float) prior a variant is associated with trait1
        p2 (float) prior a variant is associated with trait2
        p12 (float) prior a variant is associated with trait1 and trait2
    Returns:
        (array) colocalization posterior probabilities 
    """     
    lsum = l1 + l2
    lh0_abf = np.array([0])
    lh1_abf = np.log(p1) + logsumexp(l1)
    lh2_abf = np.log(p2) + logsumexp(l2)
    lh3_abf = np.log(p1) + np.log(p2) + logdiff(logsumexp(l1) + logsumexp(l2), logsumexp(lsum))
    lh4_abf = np.log(p12) + logsumexp(lsum)
    all_abf = np.concatenate([lh0_abf, lh1_abf, lh2_abf, lh3_abf, lh4_abf], axis=None)
    return np.exp(all_abf - logsumexp(all_abf))   


# other helpers
def h4_supported(h_probs: Series) -> bool:
    """ Check to see if support for H4 hypothesis, if sum of H4 is greater that sum
        of other hypotheses
    Args:
        h_probs (pandas.Series) H0-H4 hypotheses posterior probabilities
    """
    ret_val = False
    if h_probs.H4 > (h_probs.H3 + h_probs.H2 + h_probs.H1 + h_probs.H0):
        ret_val = True
    return ret_val