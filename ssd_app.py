# ssd_app.py

# This block ensures that local modules can be found reliably.
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import streamlit as st
import pandas as pd
from ssd_core import run_ssd_analysis
from database import initialize_supabase, search_chemicals_in_db, fetch_data_for_chemicals
from ui_components import create_ssd_plot
from utils import map_taxonomic_group

# --- App Configuration ---
st.set_page_config(
    page_title="SSTAC SSD Web App",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize Supabase connection
db = initialize_supabase()

# --- Initialize Session State ---
if 'data' not in st.session_state:
    st.session_state.data = None # Will hold the uploaded dataframe
if 'selected_chemicals_db' not in st.session_state:
    st.session_state.selected_chemicals_db = []
if 'selected_chemicals_csv' not in st.session_state:
    st.session_state.selected_chemicals_csv = []
if 'search_term' not in st.session_state:
    st.session_state.search_term = ""

# --- Main App UI ---
st.title("Species Sensitivity Distribution Generator")
st.markdown("A tool for ecotoxicological risk assessment using model-averaged SSDs.")

# --- Sidebar UI ---
with st.sidebar:
    st.header("Configuration")
    
    data_source = st.radio(
        "Select Data Source",
        ('Database Search', 'Upload CSV'),
        horizontal=True
    )

    # --- UI LOGIC FOR DATABASE SEARCH ---
    if data_source == 'Database Search':
        if db:
            st.session_state.search_term = st.text_input(
                "Search for a chemical", 
                value=st.session_state.search_term,
                placeholder="e.g., Copper"
            )

            search_results = []
            if st.session_state.search_term:
                search_results = search_chemicals_in_db(db, st.session_state.search_term)

            select_all = st.checkbox("Select all search results")
            
            default_selection = st.session_state.selected_chemicals_db
            if select_all:
                default_selection = search_results
            
            combined_options = sorted(list(set(search_results + st.session_state.selected_chemicals_db)))
            
            st.session_state.selected_chemicals_db = st.multiselect(
                "Select from results",
                options=combined_options,
                default=default_selection
            )

    # --- UI LOGIC FOR CSV UPLOAD ---
    else:
        uploaded_file = st.file_uploader("Upload your data", type=["csv"])
        if uploaded_file:
            st.session_state.data = pd.read_csv(uploaded_file)
            st.success("File uploaded successfully!")
        
        # --- NEW: CHEMICAL SELECTION FOR UPLOADED CSV ---
        if st.session_state.data is not None:
            if 'chemical_name' in st.session_state.data.columns:
                # Get unique chemicals from the uploaded file
                csv_chemicals = sorted(st.session_state.data['chemical_name'].unique())
                
                # Display the multiselect for the CSV
                st.session_state.selected_chemicals_csv = st.multiselect(
                    "Select Chemicals from CSV", 
                    options=csv_chemicals, 
                    default=csv_chemicals # Default to selecting all chemicals in the file
                )
            else:
                st.warning("Warning: Your CSV does not contain a 'chemical_name' column. The entire file will be analyzed as a single dataset.")
        # --- END OF NEW LOGIC ---

    st.header("SSD Parameters")
    water_type = st.radio("Media Type", ('Freshwater', 'Marine'), horizontal=True)
    agg_method = st.selectbox("Species Aggregation", ('Geometric Mean', 'Most Sensitive'))
    analysis_mode = st.selectbox("Analysis Mode", ('Model Averaging', 'Single Distribution'))

    selected_dist = None
    if analysis_mode == 'Single Distribution':
        selected_dist = st.selectbox("Select Distribution", ('Log-Normal', 'Log-Logistic', 'Weibull', 'Gamma'))

    st.header("Protection Level")
    p_value = st.slider("HCp Percentile (p-value)", 0.01, 0.50, 0.05, 0.01)
    n_boot = st.number_input("Bootstrap Iterations", 100, 5000, 1000, 100)
    
    st.header("Reproducibility")
    random_seed = st.number_input("Random Seed (for reproducible results)", 1, 999999, 42, help="Using the same seed ensures identical results across runs. Change this to get different bootstrap samples.")

    run_button = st.button("Generate SSD", type="primary")

# --- Main Panel for Results ---
if run_button:
    raw_data = None
    selected_chemicals_for_title = []

    if data_source == 'Database Search' and st.session_state.selected_chemicals_db and db:
        raw_data = fetch_data_for_chemicals(db, st.session_state.selected_chemicals_db)
        selected_chemicals_for_title = st.session_state.selected_chemicals_db

    elif data_source == 'Upload CSV' and st.session_state.data is not None:
        raw_data = st.session_state.data
        # If chemical selection is available for the CSV, filter the data
        if 'chemical_name' in raw_data.columns and st.session_state.selected_chemicals_csv:
            raw_data = raw_data[raw_data['chemical_name'].isin(st.session_state.selected_chemicals_csv)]
            selected_chemicals_for_title = st.session_state.selected_chemicals_csv
        else:
            selected_chemicals_for_title = ["Uploaded Data"]

    if raw_data is not None and not raw_data.empty:
        media_code = 'FW' if water_type == 'Freshwater' else 'MW'
        analysis_data = raw_data[raw_data['media_type'] == media_code].copy()
        
        if 'species_group' in analysis_data.columns:
            analysis_data['broad_group'] = analysis_data['species_group'].apply(map_taxonomic_group)
        
        status_box = st.status("Running SSD analysis...", expanded=True)
        progress_bar = status_box.progress(0, text="Initializing...")

        results, log_messages = run_ssd_analysis(
            data=analysis_data,
            species_col="species_scientific_name",
            value_col="conc1_mean",
            p_value=p_value,
            mode=analysis_mode.lower().replace(' ', '_'),
            selected_dist=selected_dist,
            n_boot=n_boot,
            progress_bar=progress_bar,
            random_seed=random_seed
        )
        
        status_box.update(label="Analysis complete!", state="complete", expanded=False)

        if results:
            hcp_val, ci_lower, ci_upper = results['hcp'], results['hcp_ci_lower'], results['hcp_ci_upper']
            st.success(f"**Analysis Successful:** Model-averaged HC{p_value*100:.0f} is **{hcp_val:.2f} mg/L** (95% CI: {ci_lower:.2f} - {ci_upper:.2f} mg/L)")
            
            with st.expander("View Analysis Parameters Used"):
                spec_details = {
                    "Protection Level": f"{p_value*100:.0f}% (p-value: {p_value})", "Analysis Mode": analysis_mode,
                    "Species Aggregation Method": agg_method, "Media Type Filter": water_type,
                    "Bootstrap Iterations": f"{n_boot}", "Random Seed": f"{random_seed}"
                }
                for key, value in spec_details.items(): st.markdown(f"**{key}:** `{value}`")

            tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary & Plot", "📈 Model Diagnostics", "📋 Final Data", "📝 Processing Log"])
            with tab1:
                title_chemicals = ', '.join(selected_chemicals_for_title)
                st.plotly_chart(create_ssd_plot(results['plot_data'], hcp_val, "mg/L", f"SSD for {title_chemicals} ({water_type})"), use_container_width=True)
            with tab2:
                st.dataframe(results['results_df'])
            with tab3:
                plot_df_data = results['plot_data']
                final_df = pd.DataFrame({'Species': plot_df_data['species'], 'Group': plot_df_data['groups'], 'Value (mg/L)': plot_df_data['empirical_values'], 'Percentile': plot_df_data['empirical_cdf_percent']})
                st.dataframe(final_df)
            with tab4:
                if log_messages:
                    for msg in log_messages: st.warning(msg)
                else:
                    st.info("No warnings were generated.")
        else:
            st.error("SSD analysis failed. Check the processing log.")
            if log_messages:
                for msg in log_messages: st.warning(msg)
    else:
        st.warning("Please select/upload data and click 'Generate SSD'.")