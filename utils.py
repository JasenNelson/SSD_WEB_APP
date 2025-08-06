import pandas as pd
import numpy as np

TAXONOMIC_MAPPING = {
    'Fish': ['Fish'],
    'Invertebrate': ['Invertebrate', 'Aquatic Invertebrates', 'Invertebrates', 'Crustaceans', 'Crustacean', 'Insects', 'Molluscs', 'Mollusc', 'Worms', 'Zooplankton'],
    'Plant': ['Plant', 'Algae', 'Aquatic Plants', 'Plants (Seedlings)', 'Plants', 'Algae/Plants'],
    'Amphibian': ['Amphibian', 'Amphibians']
}

def map_taxonomic_group(ecotox_group):
    """Maps detailed ECOTOX group to broader CCME category."""
    if pd.isna(ecotox_group): return "Other"
    for broad_group, detailed_list in TAXONOMIC_MAPPING.items():
        if ecotox_group in detailed_list:
            return broad_group
    return "Other"

def calculate_aicc(k, log_likelihood, n):
    """Calculates the Akaike Information Criterion corrected for small samples."""
    if n - k - 1 <= 0:
        return np.inf
    return 2 * k - 2 * log_likelihood + (2 * k**2 + 2 * k) / (n - k - 1)

def convert_df_to_csv(df):
    """Converts a DataFrame to a CSV string for downloading."""
    return df.to_csv(index=False).encode('utf-8')