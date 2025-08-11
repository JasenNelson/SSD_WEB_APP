# ssd_app.py
import streamlit as st
import pandas as pd
from ssd_core import run_ssd_analysis
from database import fetch_all_chemicals, fetch_data_for_chemical
from ui_components import build_ssd_plot
from utils import get_taxonomic_group_map

# --- App Configuration ---
st.set_page_config(
    page_title="SSTAC SSD Web App",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Initialize Session State ---
# This ensures that variables persist across user interactions.
if 'data' not in st.session_state:
    st.session_state.data = None
if 'selected_chemicals' not in st.session_state:
    st.session_state.selected_chemicals = []

# --- Main App UI ---
st.title("Species Sensitivity Distribution (SSD) Generator")
st.markdown("A tool for ecotoxicological risk assessment using model-averaged SSDs.")

# --- Sidebar UI ---
with st.sidebar:
    st.header("Configuration")
    
    data_source = st.radio(
        "Select Data Source",
        ('Database Search', 'Upload CSV'),
        horizontal=True
    )

    if data_source == 'Database Search':
        try:
            chemical_list = fetch_all_chemicals()
            st.session_state.selected_chemicals = st.multiselect(
                "Select Chemicals",
                options=chemical_list,
                default=st.session_state.selected_chemicals
            )
        except Exception as e:
            st.error(f"Failed to connect to the database: {e}")
            st.session_state.selected_chemicals = []

    else: # Upload CSV
        uploaded_file = st.file_uploader("Upload your data", type=["csv"])
        if uploaded_file:
            st.session_state.data = pd.read_csv(uploaded_file)
            st.success("File uploaded successfully!")

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

    # The main action button
    run_button = st.button("Generate SSD", type="primary")

# --- Main Panel for Results ---
# This section only runs when the button is clicked.
if run_button:
    # Clear previous results
    st.session_state.data = None
    
    # Fetch data if using database
    if data_source == 'Database Search' and st.session_state.selected_chemicals:
        try:
            st.session_state.data = fetch_data_for_chemical(st.session_state.selected_chemicals)
        except Exception as e:
            st.error(f"Failed to fetch data: {e}")

    # Check if data is available to be processed
    if st.session_state.data is not None and not st.session_state.data.empty:
        # Create a status box to show progress
        status_box = st.status("Running SSD analysis...", expanded=True)
        progress_bar = status_box.progress(0, text="Initializing...")

        # Run the core analysis function
        results, log_messages = run_ssd_analysis(
            data=st.session_state.data,
            species_col="species_scientific_name",
            value_col="conc1_mean",
            p_value=p_value,
            mode=analysis_mode.lower().replace(' ', '_'), # e.g., 'model_averaging'
            selected_dist=selected_dist,
            n_boot=n_boot,
            progress_bar=progress_bar
        )
        
        status_box.update(label="Analysis complete!", state="complete", expanded=False)

        # Display results if the analysis was successful
        if results:
            hcp_val = results['hcp']
            ci_lower = results['hcp_ci_lower']
            ci_upper = results['hcp_ci_upper']

            st.success(
                f"**Analysis Successful:** Model-averaged HC{p_value*100:.0f} is **{hcp_val:.2f} mg/L** "
                f"(95% CI: {ci_lower:.2f} - {ci_upper:.2f} mg/L)"
            )

            # --- ⬇️ BLOCK B: SPEC DISPLAY (Added) ⬇️ ---
            with st.expander("View Analysis Parameters Used"):
                spec_details = {
                    "Protection Level": f"{p_value*100:.0f}% (p-value: {p_value})",
                    "Analysis Mode": analysis_mode,
                    "Species Aggregation Method": agg_method,
                    "Media Type Filter": water_type,
                    "Bootstrap Iterations": f"{n_boot}"
                }
                for key, value in spec_details.items():
                    st.markdown(f"**{key}:** `{value}`")
            # --- ⬆️ END OF BLOCK B ⬆️ ---

            # Create tabs for organized output
            tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary & Plot", "📈 Model Diagnostics", "📋 Final Data", "📝 Processing Log"])

            with tab1:
                st.plotly_chart(build_ssd_plot(results['plot_data']), use_container_width=True)

            with tab2:
                st.dataframe(results['results_df'])
                st.download_button(
                    label="Download Diagnostics (CSV)",
                    data=results['results_df'].to_csv().encode('utf-8'),
                    file_name="gof_diagnostics.csv",
                    mime="text/csv"
                )

            with tab3:
                # Prepare final data for display and download
                plot_df_data = results['plot_data']
                final_df = pd.DataFrame({
                    'Species': plot_df_data['species'],
                    'Group': plot_df_data['groups'],
                    'Value (mg/L)': plot_df_data['empirical_values'],
                    'Percentile': plot_df_data['empirical_cdf_percent']
                })
                st.dataframe(final_df)
                st.download_button(
                    label="Download Final Data (CSV)",
                    data=final_df.to_csv().encode('utf-8'),
                    file_name="final_plot_data.csv",
                    mime="text/csv"
                )
            
            with tab4:
                if log_messages:
                    for msg in log_messages:
                        st.warning(msg)
                else:
                    st.info("No warnings were generated during processing.")
        else:
            # Display errors if analysis failed
            st.error("SSD analysis failed to produce valid results. Check the processing log for details.")
            if log_messages:
                for msg in log_messages:
                    st.warning(msg)
    else:
        st.warning("Please select chemicals or upload a file to generate an SSD.")