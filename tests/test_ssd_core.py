# tests/test_ssd_core.py

import numpy as np
import pandas as pd
import pytest
import sys
import os

# --- This helps Python find your main code files ---
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# We can now be sure these imports will work
from ssd_core import run_ssd_analysis
# --- THIS IS THE CORRECTED LINE ---
from utils import map_taxonomic_group # Use the correct function name
# --- END OF CORRECTION ---

@pytest.fixture
def boron_data():
    """A pytest fixture to load the boron dataset for tests."""
    file_path = "boron_data.csv"
    data = pd.read_csv(file_path)
    
    # The main app creates a 'broad_group' column for plotting.
    # We must replicate that step here in our test.
    # This now uses the correct mapping function.
    data['broad_group'] = data['species_group'].apply(map_taxonomic_group)
    
    return data

def test_boron_hc5_current_spec(boron_data):
    """
    Tests that the HC5 calculation for Boron is correct.
    This acts as a regression test to prevent future bugs.
    """
    np.random.seed(42) 
    
    out, logs = run_ssd_analysis(
        data=boron_data, 
        species_col="species_scientific_name", 
        value_col="conc1_mean",
        p_value=0.05, 
        mode="model_averaging", 
        n_boot=300
    )
    
    assert out is not None, "Analysis returned no results."
    assert "hcp" in out, "Results dictionary is missing 'hcp' key."
    assert 3.5 <= out["hcp"] <= 4.0, f"HC5 unexpectedly changed. Expected ~3.8, but got: {out['hcp']}"

def test_boron_hc20_sanity_check(boron_data):
    """
    Tests that the HC20 calculation matches the known benchmark (~27.5 mg/L).
    """
    np.random.seed(42)
    
    out, logs = run_ssd_analysis(
        data=boron_data, 
        species_col="species_scientific_name", 
        value_col="conc1_mean",
        p_value=0.20, 
        mode="model_averaging", 
        n_boot=300
    )

    assert out is not None, "Analysis returned no results."
    assert "hcp" in out, "Results dictionary is missing 'hcp' key."
    assert 25 <= out["hcp"] <= 29, f"HC20 unexpectedly changed. Expected ~27, but got: {out['hcp']}"
