# ssd_app.py
import streamlit as st
import pandas as pd
from database import initialize_supabase, search_chemicals_in_db, fetch_data_for_chemicals
from core_analysis import run_ssd_analysis
from ui_components import create_ssd_plot, render_diagnostics_table
from utils import convert_df_to_csv, map_taxonomic_group

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="Species Sensitivity Distribution (SSD) Generator", layout="wide")
st.title("Species Sensitivity Distribution (SSD) Generator")

# --- INITIALIZATION ---
st.session_state.setdefault('selected_chemicals', [])
supabase_client = initialize_supabase()

# --- SIDEBAR ---
with st.sidebar:
    st.header("1. Data Source")
    data_source = st.radio("Choose data source", ["Database Search", "Upload File"], horizontal=True)

    st.header("2. Chemical Selection")
    if data_source == "Database Search":
        if supabase_client:
            search_term = st.text_input("Search for a chemical in the database:", key="search_term")
            if search_term:
                st.session_state.chemical_options = search_chemicals_in_db(supabase_client, search_term)
                st.multiselect("Select chemicals for analysis:", st.session_state.chemical_options, key="selected_chemicals")
        else:
            st.warning("Database connection failed. Please check your credentials.")
    else:
        uploaded_file = st.file_uploader("Upload a CSV file", type=['csv'])
        if uploaded_file:
            df_upload = pd.read_csv(uploaded_file)
            st.session_state.chemical_options = df_upload['chemical_name'].unique().tolist()
            st.multiselect("Select chemicals from your file:", st.session_state.chemical_options, key="selected_chemicals")

    st.header("3. Guideline & Filter Options")
    water_type = st.selectbox("Filter by Water Type:", ["Freshwater (FW)", "Marine (MW)", "Both"])

    st.header("4. SSD Parameters")
    agg_method = st.selectbox("Handle multiple values per species:", ["Geometric Mean", "Most Sensitive (Minimum)"])
    analysis_mode = st.selectbox("Analysis Mode:", ["Model Averaging", "Single Distribution"])
    if analysis_mode == 'Single Distribution':
        selected_dist = st.selectbox("Select Distribution:", ['Log-Normal', 'Log-Logistic', 'Weibull', 'Gamma'])
    else:
        selected_dist = None

    st.header("5. Protection & Safety")
    hcp_percentile = st.number_input("Hazard Concentration (HCp) Percentile (%)", min_value=0.1, max_value=99.9, value=5.0, step=0.1, format="%.1f")
    n_boot = st.slider("Number of Bootstrap Iterations", min_value=100, max_value=10000, value=1000, step=100)

    generate_button = st.button("Generate SSD", type="primary", use_container_width=True)

# --- MAIN PANEL ---
if not st.session_state.selected_chemicals:
    st.info("Please select a data source and at least one chemical to begin the analysis.")
    st.stop()

if generate_button:
    if data_source == "Database Search" and supabase_client:
        df = fetch_data_for_chemicals(supabase_client, st.session_state.selected_chemicals)
    elif data_source == "Upload File" and uploaded_file is not None:
        df = df_upload.copy()
    else:
        st.error("Please select a valid data source and chemicals.")
        st.stop()
    
    if df.empty:
        st.error("No data could be retrieved for the selected chemicals. Please try a different selection.")
        st.stop()

    st.header("📈 Results")

    # CHANGED: Replaced the simple spinner with a more detailed status indicator
    with st.status("Processing data and running analysis...", expanded=True) as status:
        st.write("Filtering and preparing data...")
        proc_df = df[df['chemical_name'].isin(st.session_state.selected_chemicals)].copy()

        if 'media_type' in proc_df.columns and water_type != 'Both':
            proc_df = proc_df[proc_df['media_type'] == ('FW' if water_type == 'Freshwater (FW)' else 'MW')]
        
        # CHANGED: A more robust filter that excludes any endpoint starting with 'NR' (case-insensitive)
        if 'endpoint' in proc_df.columns:
            # The ~ symbol means "NOT", so we keep rows that DO NOT start with 'NR'.
            proc_df = proc_df[~proc_df['endpoint'].astype(str).str.upper().str.startswith('NR')]

        if proc_df.empty:
            status.update(label="Analysis Failed!", state="error", expanded=True)
            st.error("No data remains after applying filters.")
            st.stop()

        proc_df['conc1_mean'] = pd.to_numeric(proc_df['conc1_mean'], errors='coerce')
        proc_df.dropna(subset=['conc1_mean', 'species_scientific_name'], inplace=True)
        
        # Map to broad groups for plotting
        proc_df['broad_group'] = proc_df['species_group'].apply(map_taxonomic_group)
        
        # Aggregate data per species
        if agg_method == 'Geometric Mean':
            agg_func = lambda x: pd.Series({
                'conc1_mean': x['conc1_mean'].prod()**(1/len(x)),
                'endpoint': x['endpoint'].iloc[0],
                'publication_year': x['publication_year'].iloc[0],
                'author': x['author'].iloc[0],
                'title': x['title'].iloc[0],
                'chemical_name': x['chemical_name'].iloc[0],
                'broad_group': x['broad_group'].iloc[0]
            })
            final_agg_data = proc_df.groupby('species_scientific_name').apply(agg_func).reset_index()
        else: # Most Sensitive (Minimum)
            final_agg_data = proc_df.loc[proc_df.groupby('species_scientific_name')['conc1_mean'].idxmin()]
        
        st.write("Fitting distributions and starting bootstrap analysis...")
        
        # NEW: Create a progress bar to show analysis progress
        progress_bar = st.progress(0, text="Bootstrap analysis starting...")
        mode_arg = 'single' if analysis_mode == 'Single Distribution' else 'average'
        
        # CHANGED: Pass the progress_bar object to the analysis function
        results, log_messages = run_ssd_analysis(
            data=final_agg_data,
            species_col='species_scientific_name',
            value_col='conc1_mean',
            p_value=hcp_percentile / 100,
            mode=mode_arg,
            selected_dist=selected_dist,
            n_boot=n_boot,
            progress_bar=progress_bar
        )
        
        if not results:
            status.update(label="Analysis Failed!", state="error", expanded=True)
            st.error(f"SSD Calculation Error: {log_messages[0]}")
            st.stop()
        
        status.update(label="Analysis Complete!", state="complete", expanded=False)

    # --- RESULTS DISPLAY ---
    p_value = results['plot_data']['p_value']
    hcp_value = results['hcp']
    ci_lower, ci_upper = results['hcp_ci_lower'], results['hcp_ci_upper']

    tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary & Plot", "🔎 Model Diagnostics", "📋 Final Data", "⚙️ Processing Log"])

    with tab1:
        st.subheader(f"Hazard Concentration (HC{p_value*100:.0f}) Summary")
        
        col1, col2, col3 = st.columns(3)
        col1.metric(f"HC{p_value*100:.0f} Estimate", f"{hcp_value:.3g} mg/L")
        col2.metric("95% Lower CI", f"{ci_lower:.3g} mg/L")
        col3.metric("95% Upper CI", f"{ci_upper:.3g} mg/L")
        
        st.subheader("Species Sensitivity Distribution Plot")
        plot_title = f"SSD for {', '.join(st.session_state.selected_chemicals)}"
        ssd_fig = create_ssd_plot(results['plot_data'], hcp_value, "mg/L", plot_title)
        st.plotly_chart(ssd_fig, use_container_width=True)

    with tab2:
        st.subheader("Model Goodness-of-Fit Diagnostics")
        diagnostics_df = render_diagnostics_table(results['results_df'], p_value*100)
        st.download_button("Download Diagnostics (CSV)", convert_df_to_csv(diagnostics_df), "diagnostics.csv", "text/csv", use_container_width=True)

    with tab3:
        st.subheader("Aggregated Data with Source Information")
        source_cols = ['species_scientific_name', 'broad_group', 'conc1_mean', 'endpoint', 'publication_year', 'author', 'title', 'chemical_name']
        display_cols = [col for col in source_cols if col in final_agg_data.columns]
        display_df = final_agg_data[display_cols].rename(columns={'conc1_mean': 'Value (mg/L)', 'species_scientific_name': 'Species', 'broad_group': 'Group', 'chemical_name': 'Chemical'})
        st.dataframe(display_df, use_container_width=True)
        st.download_button("Download Final Data (CSV)", convert_df_to_csv(display_df), "final_data.csv", "text/csv", use_container_width=True)

    with tab4:
        st.subheader("Analysis Processing Log")
        if log_messages:
            for msg in log_messages:
                st.info(msg)
        else:
            st.success("No warnings or errors reported during analysis.")