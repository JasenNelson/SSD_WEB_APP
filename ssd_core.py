# ssd_core.py
import pandas as pd
import numpy as np
from scipy import stats
import warnings
from utils import calculate_aicc

# --- ROBUST WARNING SUPPRESSION ---
warnings.filterwarnings("ignore", category=RuntimeWarning, module='scipy')
try:
    warnings.filterwarnings("ignore", category=stats.FitConstantWarning, module='scipy')
except AttributeError:
    pass

# This structure correctly uses log-transformation for log-distributions
# and applies the necessary Jacobian correction for AICc comparison.
DISTRIBUTIONS = {
    'Log-Normal':   {'dist': stats.norm, 'k': 2, 'log': True},
    'Log-Logistic': {'dist': stats.logistic, 'k': 2, 'log': True},
    'Weibull':      {'dist': stats.weibull_min, 'k': 2, 'log': False},
    'Gamma':        {'dist': stats.gamma, 'k': 2, 'log': False},
}

def _fit_single_distribution(dist_name, model_info, data, p_value):
    """
    Fits a single distribution. Handles log-transformation and Jacobian correction
    to ensure all model AICc values are comparable.
    """
    is_log = model_info.get('log', False)
    target_data = np.log(data) if is_log else data
    
    # This ensures a correct 2-parameter fit for all models.
    params = model_info['dist'].fit(target_data, floc=0) if dist_name in ['Weibull', 'Gamma'] else model_info['dist'].fit(target_data)
    
    # Calculate log-likelihood based on the fit
    log_likelihood = np.sum(model_info['dist'].logpdf(target_data, *params))
    
    # CRITICAL: Apply Jacobian correction for models fit on log-transformed data
    if is_log:
        log_likelihood -= np.sum(np.log(data))

    aicc = calculate_aicc(model_info['k'], log_likelihood, len(data))
    
    # Calculate HCp, converting back from log-scale if necessary
    hcp = np.exp(model_info['dist'].ppf(p_value, *params)) if is_log else model_info['dist'].ppf(p_value, *params)
    
    ks_stat, ks_pvalue = stats.kstest(target_data, model_info['dist'].cdf, args=params)
    ad_stat = np.nan 

    return {'name': dist_name, 'params': params, 'aicc': aicc, 'hcp': hcp, 'ks_pvalue': ks_pvalue, 'ad_statistic': ad_stat, 'dist_obj': model_info['dist'], 'is_log': is_log}

def run_ssd_analysis(data, species_col, value_col, p_value, mode='model_averaging', selected_dist=None, n_boot=1000, progress_bar=None, random_seed=42):
    """
    Run SSD analysis with reproducible results.
    
    Args:
        data: DataFrame containing toxicity data
        species_col: Column name for species identifiers
        value_col: Column name for concentration values
        p_value: Protection level (e.g., 0.05 for HC5)
        mode: Analysis mode ('model_averaging' or 'single_distribution')
        selected_dist: Specific distribution for single mode
        n_boot: Number of bootstrap iterations
        progress_bar: Streamlit progress bar for UI feedback
        random_seed: Random seed for reproducible bootstrap results (default: 42)
    """
    # Set random seed for reproducible results - this ensures the same bootstrap
    # samples are generated every time, leading to consistent confidence intervals
    np.random.seed(random_seed)
    
    valid_data_df = data[data[value_col] > 0].copy()
    valid_data = valid_data_df[value_col]
    if len(valid_data) < 5: return None, ["Not enough valid data points (minimum 5 required)."]
    
    n = len(valid_data); log_messages = []
    
    dists_to_fit = [selected_dist] if mode == 'single_distribution' else list(DISTRIBUTIONS.keys())
    model_fits = []
    for name in dists_to_fit:
        if name not in DISTRIBUTIONS: continue
        try:
            fit = _fit_single_distribution(name, DISTRIBUTIONS[name], valid_data, p_value)
            if not fit or not np.isfinite(fit['hcp']) or fit['hcp'] <= 0:
                log_messages.append(f"Warning: Could not derive a valid positive HCp for the '{name}' distribution. It will be excluded.")
                continue
            model_fits.append(fit)
        except Exception as e:
            log_messages.append(f"Warning: Could not fit '{name}' to the original data. It will be excluded. Details: {e}")
    if not model_fits: return None, ["Failed to fit any valid distributions to the data."]
    
    results_df = pd.DataFrame(model_fits)
    
    # --- DEBUGGING CANARY ---
    # This print statement will ONLY appear in the logs if the new code is running.
    print("--- DEBUG v3: AICc values BEFORE weighting ---")
    print(results_df[['name', 'aicc']])
    print("--------------------------------------------")
    # --- END DEBUGGING ---
    
    results_df['weight'] = 1.0 if mode == 'single_distribution' else np.exp(-0.5 * (results_df['aicc'] - results_df['aicc'].min())) / np.sum(np.exp(-0.5 * (results_df['aicc'] - results_df['aicc'].min())))
    final_hcp = np.sum(results_df['weight'] * results_df['hcp'])

    # Standard non-parametric bootstrap
    boot_hcps, boot_cdfs = [], []
    # Extend x-range to ensure it covers the HC5 area and beyond
    min_val = valid_data.min()
    max_val = valid_data.max()
    x_range_log = np.linspace(np.log(min_val) - 2, np.log(max_val) + 2, 200)
    
    for i in range(n_boot):
        try:
            boot_sample = valid_data.sample(n=n, replace=True)
            
            b_fits = []
            for name in dists_to_fit:
                if name not in DISTRIBUTIONS: continue
                try:
                    b_fit = _fit_single_distribution(name, DISTRIBUTIONS[name], boot_sample, p_value)
                    if b_fit and np.isfinite(b_fit['hcp']) and b_fit['hcp'] > 0:
                        b_fits.append(b_fit)
                except Exception: continue
            
            if not b_fits: continue
            
            b_results_df = pd.DataFrame(b_fits)
            b_results_df['weight'] = 1.0 if mode == 'single_distribution' else np.exp(-0.5 * (b_results_df['aicc'] - b_results_df['aicc'].min())) / np.sum(np.exp(-0.5 * (b_results_df['aicc'] - b_results_df['aicc'].min())))
            boot_hcps.append(np.sum(b_results_df['weight'] * b_results_df['hcp']))
            
            b_avg_cdf = np.zeros_like(x_range_log)
            for _, row in b_results_df.iterrows():
                if row['is_log']:
                    # For log-transformed distributions, use log-transformed x-range
                    cdf_vals = row['dist_obj'].cdf(x_range_log, *row['params'])
                else:
                    # For non-log distributions, use original scale x-range
                    cdf_vals = row['dist_obj'].cdf(np.exp(x_range_log), *row['params'])
                b_avg_cdf += row['weight'] * cdf_vals
            boot_cdfs.append(b_avg_cdf)
            
            if progress_bar:
                progress_bar.progress((i + 1) / n_boot, text=f"Running bootstrap iteration {i+1} of {n_boot}")
        except Exception as e:
            log_messages.append(f"Warning: A bootstrap iteration failed unexpectedly. Details: {e}")
            continue

    if len(boot_hcps) < n_boot * 0.8:
        log_messages.append(f"Warning: Confidence intervals may be unreliable. Only {len(boot_hcps)}/{n_boot} successful bootstrap iterations.")
    
    hcp_ci_lower, hcp_ci_upper = np.percentile(boot_hcps, [2.5, 97.5]) if len(boot_hcps) > 2 else (np.nan, np.nan)
    lower_ci_curve, upper_ci_curve = np.percentile(boot_cdfs, [2.5, 97.5], axis=0) if len(boot_hcps) > 2 else (np.full_like(x_range_log, np.nan), np.full_like(x_range_log, np.nan))
    
    final_avg_cdf = np.zeros_like(x_range_log)
    for _, row in results_df.iterrows():
        if row['is_log']:
            # For log-transformed distributions, use log-transformed x-range
            cdf_vals = row['dist_obj'].cdf(x_range_log, *row['params'])
        else:
            # For non-log distributions, use original scale x-range
            cdf_vals = row['dist_obj'].cdf(np.exp(x_range_log), *row['params'])
        final_avg_cdf += row['weight'] * cdf_vals

    full_data = data.sort_values(by=value_col)
    empirical_cdf_percent = (np.arange(1, len(full_data) + 1) / (len(full_data) + 1)) * 100
    
    # Cap the fitted curve at the empirical maximum to prevent overshooting
    empirical_max = np.max(empirical_cdf_percent)
    final_avg_cdf = np.minimum(final_avg_cdf * 100, empirical_max)
    lower_ci_curve = np.minimum(lower_ci_curve * 100, empirical_max)
    upper_ci_curve = np.minimum(upper_ci_curve * 100, empirical_max)
    
    plot_data = {'empirical_values': full_data[value_col].values, 'empirical_cdf_percent': empirical_cdf_percent, 'fitted_x_range': np.exp(x_range_log), 'fitted_y_cdf_percent': final_avg_cdf, 'lower_ci_percent': lower_ci_curve, 'upper_ci_percent': upper_ci_curve, 'species': full_data[species_col].tolist(), 'groups': full_data['broad_group'].tolist(), 'p_value': p_value}
    final_results = {'hcp': final_hcp, 'hcp_ci_lower': hcp_ci_lower, 'hcp_ci_upper': hcp_ci_upper, 'results_df': results_df, 'plot_data': plot_data}
    return final_results, log_messages