# core_analysis.py
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

# --- DEFINITIVE FIX: Use SciPy's correct name for the Log-Logistic distribution ('fisk') ---
DISTRIBUTIONS = {
    'Log-Normal':   {'dist': stats.lognorm, 'k': 2},
    'Log-Logistic': {'dist': stats.fisk, 'k': 2},
    'Weibull':      {'dist': stats.weibull_min, 'k': 2},
    'Gamma':        {'dist': stats.gamma, 'k': 2},
}

def _fit_single_distribution(dist_name, model_info, data, p_value):
    """Fits a single distribution directly on the original concentration data."""
    dist_obj = model_info['dist']
    
    # Fit the distribution. Use floc=0 for distributions that require it.
    params = dist_obj.fit(data, floc=0) if dist_name in ['Weibull', 'Gamma'] else dist_obj.fit(data)
    
    # Calculate log-likelihood on the original data scale for all models.
    log_likelihood = np.sum(dist_obj.logpdf(data, *params))
    
    # AICc values are now directly comparable.
    aicc = calculate_aicc(model_info['k'], log_likelihood, len(data))
    
    # Calculate HCp directly.
    hcp = dist_obj.ppf(p_value, *params)
    
    # Run goodness-of-fit tests.
    ks_stat, ks_pvalue = stats.kstest(data, lambda x: dist_obj.cdf(x, *params))
    
    ad_stat = np.nan 

    return {'name': dist_name, 'params': params, 'aicc': aicc, 'hcp': hcp, 'ks_pvalue': ks_pvalue, 'ad_statistic': ad_stat, 'dist_obj': dist_obj}

def run_ssd_analysis(data, species_col, value_col, p_value, mode='average', selected_dist=None, n_boot=1000, progress_bar=None):
    valid_data_df = data[data[value_col] > 0].copy()
    valid_data = valid_data_df[value_col]
    if len(valid_data) < 5: return None, ["Not enough valid data points (minimum 5 required)."]
    
    n = len(valid_data); log_messages = []
    
    dists_to_fit = [selected_dist] if mode == 'single' else list(DISTRIBUTIONS.keys())
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
    results_df['weight'] = 1.0 if mode == 'single' else np.exp(-0.5 * (results_df['aicc'] - results_df['aicc'].min())) / np.sum(np.exp(-0.5 * (results_df['aicc'] - results_df['aicc'].min())))
    final_hcp = np.sum(results_df['weight'] * results_df['hcp'])

    # --- SIMPLIFIED & ROBUST BOOTSTRAP ---
    boot_hcps, boot_cdfs = [], []
    x_range = np.logspace(np.log10(valid_data.min() * 0.5), np.log10(valid_data.max() * 1.2), 200)
    
    for i in range(n_boot):
        try:
            chosen_model_row = results_df.sample(n=1, weights='weight').iloc[0]
            
            # 1. Generate new data on the original scale.
            boot_sample = chosen_model_row['dist_obj'].rvs(*chosen_model_row['params'], size=n)
            boot_sample = boot_sample[np.isfinite(boot_sample) & (boot_sample > 0)]
            if len(boot_sample) < 5: continue

            # 2. Re-fit model(s) on the clean bootstrap sample
            b_fits = []
            for name in dists_to_fit:
                if name not in DISTRIBUTIONS: continue
                try:
                    b_fit = _fit_single_distribution(name, DISTRIBUTIONS[name], boot_sample, p_value)
                    if b_fit and np.isfinite(b_fit['hcp']) and b_fit['hcp'] > 0:
                        b_fits.append(b_fit)
                except Exception: continue
            
            if not b_fits: continue
            
            # 3. Calculate and append results
            b_results_df = pd.DataFrame(b_fits)
            b_results_df['weight'] = 1.0 if mode == 'single' else np.exp(-0.5 * (b_results_df['aicc'] - b_results_df['aicc'].min())) / np.sum(np.exp(-0.5 * (b_results_df['aicc'] - b_results_df['aicc'].min())))
            boot_hcps.append(np.sum(b_results_df['weight'] * b_results_df['hcp']))
            
            b_avg_cdf = np.zeros_like(x_range)
            for _, row in b_results_df.iterrows():
                b_avg_cdf += row['weight'] * row['dist_obj'].cdf(x_range, *row['params'])
            boot_cdfs.append(b_avg_cdf)
            
            if progress_bar:
                progress_bar.progress((i + 1) / n_boot, text=f"Running bootstrap iteration {i+1} of {n_boot}")
        except Exception as e:
            log_messages.append(f"Warning: A bootstrap iteration failed unexpectedly. Details: {e}")
            continue

    if len(boot_hcps) < n_boot * 0.1: # Check if a significant number of bootstraps failed
        log_messages.append(f"Warning: Confidence intervals may be unreliable. Only {len(boot_hcps)}/{n_boot} bootstrap iterations succeeded.")
    
    hcp_ci_lower, hcp_ci_upper = np.percentile(boot_hcps, [2.5, 97.5]) if len(boot_hcps) > 2 else (np.nan, np.nan)
    lower_ci_curve, upper_ci_curve = np.percentile(boot_cdfs, [2.5, 97.5], axis=0) if len(boot_hcps) > 2 else (np.full_like(x_range, np.nan), np.full_like(x_range, np.nan))
    
    final_avg_cdf = np.zeros_like(x_range)
    for _, row in results_df.iterrows():
        final_avg_cdf += row['weight'] * row['dist_obj'].cdf(x_range, *row['params'])

    full_data = data.sort_values(by=value_col)
    plot_data = {'empirical_values': full_data[value_col].values, 'empirical_cdf_percent': (np.arange(1, len(full_data) + 1) / (len(full_data) + 1)) * 100, 'fitted_x_range': x_range, 'fitted_y_cdf_percent': final_avg_cdf * 100, 'lower_ci_percent': lower_ci_curve * 100, 'upper_ci_percent': upper_ci_curve * 100, 'species': full_data[species_col].tolist(), 'groups': full_data['broad_group'].tolist(), 'p_value': p_value}
    final_results = {'hcp': final_hcp, 'hcp_ci_lower': hcp_ci_lower, 'hcp_ci_upper': hcp_ci_upper, 'results_df': results_df, 'plot_data': plot_data}
    return final_results, log_messages