# ssd_app.py
import streamlit as st
import pandas as pd
from scipy import stats
import numpy as np

from utils import map_taxonomic_group, convert_df_to_csv
from database import initialize_supabase, search_chemicals_in_db, fetch_data_for_chemicals
from core_analysis import run_ssd_analysis, DISTRIBUTIONS
from ui_components import create_ssd_plot, render_diagnostics_table

st.set_page_config(layout="wide", page_title="Professional SSD Generator")

if 'search_results' not in st.session_state: st.session_state.search_results = []
if 'selected_chemicals' not in st.session_state: st.session_state.selected_chemicals = []

st.title("Professional Species Sensitivity Distribution (SSD) Generator")
st.markdown("A `ssdtools`-aligned analysis tool with model averaging and power-user features.")

supabase_conn = initialize_supabase()

with st.sidebar:
    st.header("⚙️ Configuration")
    data_source = st.radio("1. Data Source", ("Database", "File Upload"), horizontal=True, key="data_source")
    
    st.subheader("2. Chemical Selection")
    if data_source == "Database":
        if supabase_conn:
            search_term = st.text_input("Search for chemicals")
            if st.button("Search"):
                with st.spinner("Searching..."):
                    st.session_state.search_results = search_chemicals_in_db(supabase_conn, search_term)
            
            options = sorted(list(set(st.session_state.selected_chemicals + st.session_state.search_results)))
            
            st.write("Select chemicals from the list below:")
            select_all = st.checkbox("Select all search results")
            if select_all:
                st.session_state.selected_chemicals = st.multiselect("Selected Chemicals", options=options, default=options)
            else:
                st.session_state.selected_chemicals = st.multiselect("Selected Chemicals", options=options, default=st.session_state.selected_chemicals)
        else: st.error("Database connection unavailable.")
    else: # File Upload
        uploaded_file = st.file_uploader("Upload Data (CSV)", type=["csv"])
        if uploaded_file:
            try:
                temp_df = pd.read_csv(uploaded_file, encoding='utf-8', encoding_errors='replace')
                if 'chemical_name' in temp_df.columns:
                    chem_options = sorted(temp_df['chemical_name'].dropna().unique())
                    select_all_file = st.checkbox("Select all chemicals from file")
                    if select_all_file:
                        st.session_state.selected_chemicals = st.multiselect("Selected Chemicals", options=chem_options, default=chem_options)
                    else:
                        st.session_state.selected_chemicals = st.multiselect("Selected Chemicals", options=chem_options, default=st.session_state.selected_chemicals)
                else: st.warning("File must contain a 'chemical_name' column.")
            except Exception as e: st.error(f"Error reading file: {e}")

    st.subheader("3. Analysis Options")
    analysis_mode = st.radio("Analysis Mode", ('Model Averaging (Recommended)', 'Single Distribution'), horizontal=True)
    selected_dist = st.selectbox("Select Distribution", options=list(DISTRIBUTIONS.keys())) if analysis_mode == 'Single Distribution' else None
    
    with st.expander("Advanced Analysis Settings"):
        n_boot = st.slider("Bootstrap Iterations", min_value=100, max_value=5000, value=1000, step=100, help="Fewer iterations are faster but less precise.")
        hcp_percentile = st.number_input("HCp Percentile (%)", 1.0, 50.0, 5.0, 0.5, "%.1f")
        water_type = st.radio("Water Type", ('Freshwater (FW)', 'Marine Water (MW)', 'Both'), horizontal=True)
        data_handling = st.radio("Handle Multiple Values per Species", ('Use Geometric Mean', 'Use Most Sensitive'), horizontal=True, index=0)
    
    generate_button = st.button("Generate SSD", type="primary", use_container_width=True)

if generate_button:
    if not st.session_state.selected_chemicals: st.warning("Please select at least one chemical."); st.stop()
    
    df = None
    if st.session_state.data_source == "File Upload":
        if 'uploaded_file' in locals() and uploaded_file:
            uploaded_file.seek(0); df = pd.read_csv(uploaded_file, encoding='utf-8', encoding_errors='replace')
        else: st.warning("Please upload a file."); st.stop()
    else:
        df = fetch_data_for_chemicals(supabase_conn, st.session_state.selected_chemicals)
    if df.empty: st.error("No data could be loaded for the selected chemical(s)."); st.stop()

    with st.spinner("Processing data and running analysis..."):
        proc_df = df[df['chemical_name'].isin(st.session_state.selected_chemicals)].copy()
        if 'media_type' in proc_df.columns and water_type != 'Both':
            proc_df = proc_df[proc_df['media_type'] == ('FW' if water_type == 'Freshwater (FW)' else 'MW')]
        if proc_df.empty: st.error("No data remains after applying filters."); st.stop()
        proc_df['conc1_mean'] = pd.to_numeric(proc_df['conc1_mean'], errors='coerce')
        proc_df.dropna(subset=['conc1_mean', 'species_scientific_name'], inplace=True)
        if 'publication_year' in proc_df.columns:
            proc_df['publication_year'] = pd.to_numeric(proc_df['publication_year'], errors='coerce')
        else:
            proc_df['publication_year'] = np.nan
        proc_df['broad_group'] = proc_df['species_group'].apply(map_taxonomic_group)
        if data_handling == 'Use Geometric Mean':
            proc_df['log_conc'] = np.log(proc_df['conc1_mean'])
            gmeans = proc_df.groupby('species_scientific_name')['log_conc'].mean().apply(np.exp)
            latest_idx = proc_df.groupby('species_scientific_name')['publication_year'].idxmax().dropna()
            source_info = proc_df.loc[latest_idx].set_index('species_scientific_name')
            final_agg_data = source_info.copy()
            final_agg_data['conc1_mean'] = gmeans
            final_agg_data.reset_index(inplace=True)
        else:
            latest_idx = proc_df.groupby('species_scientific_name')['conc1_mean'].idxmin().dropna()
            final_agg_data = proc_df.loc[latest_idx]
        final_agg_data.dropna(subset=['conc1_mean'], inplace=True)

        mode_arg = 'single' if analysis_mode == 'Single Distribution' else 'average'
        results, log_messages = run_ssd_analysis(data=final_agg_data, species_col='species_scientific_name', value_col='conc1_mean', p_value=hcp_percentile / 100, mode=mode_arg, selected_dist=selected_dist, n_boot=n_boot)
        if not results: st.error(f"SSD Calculation Error: {log_messages[0]}"); st.stop()

    st.header("📈 Results")
    chemical_str = ', '.join(st.session_state.selected_chemicals)
    plot_title = f"SSD for {chemical_str} ({analysis_mode.split(' ')[0]})"
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary & Plot", "🔎 Model Diagnostics", "📋 Final Data", "⚙️ Processing Log"])
    with tab1:
        st.subheader("Hazard Concentration Summary")
        col1, col2 = st.columns(2)
        col1.metric(label=f"HC{hcp_percentile:.1f}", value=f"{results['hcp']:.4g} mg/L")
        col2.metric(label=f"95% Confidence Interval", value=f"{results['hcp_ci_lower']:.4g} – {results['hcp_ci_upper']:.4g} mg/L")
        
        st.subheader("Species Sensitivity Distribution Plot")
        ssd_fig = create_ssd_plot(results['plot_data'], results['hcp'], 'mg/L', plot_title)
        st.plotly_chart(ssd_fig, use_container_width=True)
        
        # --- RESILIENT IMAGE EXPORT ---
        try:
            img_bytes = ssd_fig.to_image(format="png", width=1200, height=700, scale=2)
            st.download_button(
                label="📥 Download Plot (PNG)", 
                data=img_bytes, 
                file_name="ssd_plot.png", 
                mime="image/png"
            )
        except Exception as e:
            st.warning(
                "Could not generate plot for download. This may be due to a temporary environment issue. "
                "Please try rebooting the app from the Streamlit Cloud dashboard menu (☰).", 
                icon="⚠️"
            )
            st.caption(f"Details: {e}")
    with tab2:
        # --- BUG FIX: Changed results['df'] to results['results_df'] ---
        diagnostics_df = render_diagnostics_table(results['results_df'], hcp_percentile)
        csv_data = convert_df_to_csv(diagnostics_df)
        st.download_button("📥 Export Diagnostics to CSV", data=csv_data, file_name="ssd_diagnostics.csv", mime="text/csv")
    with tab3:
        st.subheader("Aggregated Data with Source Information")
        source_cols = ['species_scientific_name', 'broad_group', 'conc1_mean', 'endpoint', 'publication_year', 'author', 'title', 'chemical_name']
        display_cols = [col for col in source_cols if col in final_agg_data.columns]
        display_df = final_agg_data[display_cols].rename(columns={'conc1_mean': 'Value (mg/L)', 'species_scientific_name': 'Species', 'broad_group': 'Group', 'chemical_name': 'Chemical'})
        st.dataframe(display_df, use_container_width=True)
        csv_data = convert_df_to_csv(display_df)
        st.download_button("📥 Export Final Data to CSV", data=csv_data, file_name="final_data.csv", mime="text/csv")
    with tab4:
        st.subheader("Analysis Log")
        if log_messages:
            for msg in log_messages: st.warning(msg)
        else: st.success("Analysis completed successfully.")