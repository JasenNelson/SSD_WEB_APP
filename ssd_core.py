# ssd_core.py
import pandas as pd
import numpy as np
from scipy import stats
import warnings
from utils import calculate_aicc

warnings.filterwarnings("ignore", category=RuntimeWarning, module='scipy')
try:
    warnings.filterwarnings("ignore", category=stats.FitConstantWarning, module='scipy')
except AttributeError:
    pass

DISTRIBUTIONS = {
    'Log-Normal':   {'dist': stats.norm, 'k': 2, 'log': True},
    'Log-Logistic': {'dist': stats.logistic, 'k': 2, 'log': True},
    'Weibull':      {'dist': stats.weibull_min, 'k': 2, 'log': False},
    'Gamma':        {'dist': stats.gamma, 'k': 2, 'log': False},
}

def _fit_single_distribution(dist_name, model_info, data, p_value):
    """
    Fits a single statistical distribution to the data, calculates goodness-of-fit metrics,
    and estimates the Hazard Concentration (HCp).
    """
    is_log = model_info.get('log', False)
    target_data = np.log(data) if is_log else data
    dist_obj = model_info['dist']
    
    # --- BUG FIX START ---
    # The original generic parameter handling was replaced with an explicit, robust
    # approach for each distribution type to prevent incorrect parameter passing.
    
    hcp = np.nan
    log_likelihood = np.nan
    params = ()

    if dist_name in ['Log-Normal', 'Log-Logistic']:
        # These models are fitted to log-transformed data. They have 2 parameters (loc, scale).
        loc, scale = dist_obj.fit(target_data)
        params = (loc, scale)
        log_likelihood = np.sum(dist_obj.logpdf(target_data, loc=loc, scale=scale))
        # Apply Jacobian correction for log-transformation to ensure AICc is comparable
        log_likelihood -= np.sum(np.log(data))
        hcp = np.exp(dist_obj.ppf(p_value, loc=loc, scale=scale))

    elif dist_name in ['Weibull', 'Gamma']:
        # These models are fitted to the original data as 2-parameter distributions (shape, scale)
        # by fixing the location parameter (loc) to 0.
        shape, loc, scale = dist_obj.fit(target_data, floc=0)
        params = (shape, loc, scale)
        log_likelihood = np.sum(dist_obj.logpdf(target_data, shape, loc=loc, scale=scale))
        hcp = dist_obj.ppf(p_value, shape, loc=loc, scale=scale)

    else:
        raise ValueError(f"Distribution '{dist_name}' not recognized for fitting.")

    # --- BUG FIX END ---
        
    aicc = calculate_aicc(model_info['k'], log_likelihood, len(data))
    
    # Perform goodness-of-fit tests
    ks_stat, ks_pvalue = stats.kstest(target_data, dist_obj.cdf, args=params)
    ad_stat = np.nan # Anderson-Darling requires more specific implementation per distribution

    return {
        'name': dist_name, 'params': params, 'aicc': aicc, 'hcp': hcp, 
        'ks_pvalue': ks_pvalue, 'ad_statistic': ad_stat, 
        'dist_obj': dist_obj, 'is_log': is_log
    }

def run_ssd_analysis(data, species_col, value_col, p_value, mode='model_averaging', selected_dist=None, n_boot=1000, progress_bar=None, random_seed=42):
    """
    Main function to run the Species Sensitivity Distribution analysis.
    """
    np.random.seed(random_seed)
    valid_data_df = data[data[value_col] > 0].copy()
    valid_data = valid_data_df[value_col]
    
    if len(valid_data) < 5:
        return None, ["Not enough valid data points (minimum 5 required)."]
    
    n = len(valid_data)
    log_messages = []
    
    # Determine which distributions to fit
    dists_to_fit = [selected_dist] if mode == 'single_distribution' and selected_dist else list(DISTRIBUTIONS.keys())
    
    model_fits = []
    for i, name in enumerate(dists_to_fit):
        if name not in DISTRIBUTIONS:
            continue
        try:
            if progress_bar:
                progress_bar.progress((i + 1) / len(dists_to_fit), text=f"Fitting {name} distribution...")
            
            fit = _fit_single_distribution(name, DISTRIBUTIONS[name], valid_data, p_value)
            
            if not fit or not np.isfinite(fit['hcp']) or fit['hcp'] <= 0:
                log_messages.append(f"Warning: Could not derive a valid HCp for '{name}'. Excluding.")
                continue
            
            model_fits.append(fit)
        except Exception as e:
            log_messages.append(f"Warning: Could not fit '{name}'. Excluding. Details: {e}")
            
    if not model_fits:
        return None, ["Failed to fit any valid distributions."]
    
    results_df = pd.DataFrame(model_fits).sort_values(by='aicc').reset_index(drop=True)
    
    if mode == 'model_averaging':
        # Calculate Akaike weights
        delta_aicc = results_df['aicc'] - results_df['aicc'].min()
        exp_delta = np.exp(-0.5 * delta_aicc)
        results_df['weight'] = exp_delta / np.sum(exp_delta)
        final_hcp = np.sum(results_df['weight'] * results_df['hcp'])
    else: # Single distribution mode
        results_df['weight'] = 1.0
        final_hcp = results_df['hcp'][0]

    # [Bootstrap logic would follow here...]
    # For now, returning dummy CI values as in the original snippet.
    
    return {
        'hcp': final_hcp, 
        'hcp_ci_lower': 0, 
        'hcp_ci_upper': 0, 
        'results_df': results_df, 
        'plot_data': None
    }, log_messages