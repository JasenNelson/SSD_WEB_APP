# core_analysis.py
import pandas as pd
import numpy as np
from scipy import stats
import streamlit as st
import warnings
from utils import calculate_aicc

# --- ROBUST WARNING SUPPRESSION ---
warnings.filterwarnings("ignore", category=RuntimeWarning, module='scipy')
try:
    warnings.filterwarnings("ignore", category=stats.FitConstantWarning, module='scipy')
except AttributeError:
    pass

DISTRIBUTIONS = {
    'Log-Normal':   {'dist': stats.norm, 'k': 2, 'log': True, 'scipy_name': 'norm'},
    'Log-Logistic': {'dist': stats.logistic, 'k': 2, 'log': True, 'scipy_name': 'logistic'},
    'Weibull':      {'dist': stats.weibull_min, 'k': 2, 'log': False, 'scipy_name': 'weibull_min'},
    'Gamma':        {'dist': stats.gamma, 'k': 2, 'log': False, 'scipy_name': 'gamma'},
}

def _fit_single_distribution(dist_name, model_info, data, p_value):
    """Fits a single distribution. Expects data on the original concentration scale."""
    target_data = np.log(data) if model_info['log'] else data
    params = model_info['dist'].fit(target_data, floc=0) if dist_name == 'Weibull' else model_info['dist'].fit(target_data)
    log_likelihood = np.sum(model_info['dist'].logpdf(target_data, *params))
    aicc = calculate_aicc(model_info['k'], log_likelihood, len(target_data))
    hcp = np.exp(model_info['dist'].ppf(p_value, *params)) if model_info['log'] else model_info['dist'].ppf(p_value, *params)
    ks_stat, ks_pvalue = stats.kstest(target_data, model_info['dist'].cdf, args=params)
    supported_ad = ['norm', 'logistic', 'weibull_min']
    ad_stat = stats.anderson(target_data, dist=model_info['scipy_name'])[0] if model_info['scipy_name'] in supported_ad else np.nan
    return {'name': dist_name, 'params': params, 'aicc': aicc, 'hcp': hcp, 'ks_pvalue': ks_pvalue, 'ad_statistic': ad_stat, 'dist_obj': model_info['dist'], 'is_log': model_info['log']}

@st.cache_data(show_spinner=False)
def run_ssd_analysis(data, species_col, value_col, p_value, mode='average', selected_dist=None, n_boot=1000):
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
            model_fits.append(fit)
        except Exception as e:
            log_messages.append(f"Warning: Could not fit '{name}' to the original data. It will be excluded. Details: {e}")
    if not model_fits: return None, ["Failed to fit any distributions to the original data."]
    
    results_df = pd.DataFrame(model_fits)
    results_df['weight'] = 1.0 if mode == 'single' else np.exp(-0.5 * (results_df['aicc'] - results_df['aicc'].min())) / np.sum(np.exp(-0.5 * (results_df['aicc'] - results_df['aicc'].min())))
    final_hcp = np.sum(results_df['weight'] * results_df['hcp'])

    # --- DEFINITIVE FIX for BOOTSTRAPPING LOGIC ---
    boot_hcps, boot_cdfs = [], []
    x_range_log = np.linspace(np.log(valid_data.min()) * 0.9, np.log(valid_data.max()) * 1.1, 200)
    
    for i in range(n_boot):
        try:
            chosen_model_row = results_df.sample(n=1, weights='weight').iloc[0]
            
            # 1. Generate new data.
            boot_sample = chosen_model_row['dist_obj'].rvs(*chosen_model_row['params'], size=n)
            
            # 2. Convert to original scale and apply robust cleaning.
            if chosen_model_row['is_log']:
                boot_sample = np.exp(boot_sample)
            
            boot_sample = boot_sample[np.isfinite(boot_sample)] # Remove inf/nan
            boot_sample = boot_sample[boot_sample > 0]           # Remove non-positives
            
            if len(boot_sample) < 5: continue # Skip if sample is too small after cleaning

            # 3. Re-fit model(s) on the clean bootstrap sample
            b_fits = []
            for name in dists_to_fit:
                if name not in DISTRIBUTIONS: continue
                try:
                    b_fit = _fit_single_distribution(name, DISTRIBUTIONS[name], boot_sample, p_value)
                    b_fits.append(b_fit)
                except Exception:
                    continue # Ignore if a single dist fails to fit to the bootstrap sample
            
            if len(b_fits) < len(dists_to_fit):
                log_messages.append(f"Note: One or more distributions failed to fit during bootstrap iteration {i+1}.")
            if not b_fits: continue
            
            # 4. Calculate and append results for this successful iteration
            b_results_df = pd.DataFrame(b_fits)
            b_results_df['weight'] = 1.0 if mode == 'single' else np.exp(-0.5 * (b_results_df['aicc'] - b_results_df['aicc'].min())) / np.sum(np.exp(-0.5 * (b_results_df['aicc'] - b_results_df['aicc'].min())))
            boot_hcps.append(np.sum(b_results_df['weight'] * b_results_df['hcp']))
            
            b_avg_cdf = np.zeros_like(x_range_log)
            for _, row in b_results_df.iterrows():
                cdf_vals = row['dist_obj'].cdf(x_range_log, *row['params']) if row['is_log'] else row['dist_obj'].cdf(np.exp(x_range_log), *row['params'])
                b_avg_cdf += row['weight'] * cdf_vals
            boot_cdfs.append(b_avg_cdf)
        except Exception as e:
            # Catch any unexpected errors from the main loop
            log_messages.append(f"Warning: A bootstrap iteration failed unexpectedly. Details: {e}")
            continue

    if len(boot_hcps) < 2:
        log_messages.append("Warning: Confidence intervals could not be calculated (fewer than 2 successful bootstrap iterations).")
        hcp_ci_lower, hcp_ci_upper = np.nan, np.nan
        lower_ci_curve, upper_ci_curve = np.full_like(x_range_log, np.nan), np.full_like(x_range_log, np.nan)
    else:
        hcp_ci_lower, hcp_ci_upper = np.percentile(boot_hcps, 2.5), np.percentile(boot_hcps, 97.5)
        lower_ci_curve, upper_ci_curve = np.percentile(boot_cdfs, 2.5, axis=0), np.percentile(boot_cdfs, 97.5, axis=0)
    
    final_avg_cdf = np.zeros_like(x_range_log)
    for _, row in results_df.iterrows():
        cdf_vals = row['dist_obj'].cdf(x_range_log, *row['params']) if row['is_log'] else row['dist_obj'].cdf(np.exp(x_range_log), *row['params'])
        final_avg_cdf += row['weight'] * cdf_vals

    full_data = data.sort_values(by=value_col)
    plot_data = {'empirical_values': full_data[value_col].values, 'empirical_cdf_percent': (np.arange(1, len(full_data) + 1) / (len(full_data) + 1)) * 100, 'fitted_x_range': np.exp(x_range_log), 'fitted_y_cdf_percent': final_avg_cdf * 100, 'lower_ci_percent': lower_ci_curve * 100, 'upper_ci_percent': upper_ci_curve * 100, 'species': full_data[species_col].tolist(), 'groups': full_data['broad_group'].tolist(), 'p_value': p_value}
    final_results = {'hcp': final_hcp, 'hcp_ci_lower': hcp_ci_lower, 'hcp_ci_upper': hcp_ci_upper, 'results_df': results_df, 'plot_data': plot_data}
    return final_results, log_messages