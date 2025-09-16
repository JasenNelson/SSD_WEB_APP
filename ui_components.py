# ui_components.py
import plotly.graph_objects as go
import pandas as pd
import streamlit as st

try:
    import matplotlib
    MATPLOTLIB_INSTALLED = True
except ImportError:
    MATPLOTLIB_INSTALLED = False

MARKER_STYLES = {'Fish': {'symbol': 'circle', 'color': '#1f77b4'}, 'Invertebrate': {'symbol': 'diamond', 'color': '#ff7f0e'}, 'Plant': {'symbol': 'square', 'color': '#2ca02c'}, 'Amphibian': {'symbol': 'cross', 'color': '#d62728'}, 'Other': {'symbol': 'x', 'color': '#9467bd'}}

def create_ssd_plot(plot_data, hcp, unit, title):
    if plot_data is None: return go.Figure()
    fig = go.Figure()
    
    empirical_x, empirical_y = plot_data['empirical_values'], plot_data['empirical_cdf_percent']
    species_names, species_groups = plot_data['species'], plot_data['groups']
    fitted_x, fitted_y = plot_data['fitted_x_range'], plot_data['fitted_y_cdf_percent']
    lower_ci, upper_ci = plot_data['lower_ci_percent'], plot_data['upper_ci_percent']
    
    # --- The Definitive Method for Plotting a Shaded Confidence Band ---
    fig.add_trace(go.Scatter(
        x=list(fitted_x) + list(fitted_x[::-1]),
        y=list(upper_ci) + list(lower_ci[::-1]),
        fill='toself',
        fillcolor='rgba(0,100,80,0.2)',
        line=dict(color='rgba(255,255,255,0)'),
        hoverinfo="none",
        showlegend=True,
        name='95% Confidence Interval'
    ))
    
    fig.add_trace(go.Scatter(x=fitted_x, y=fitted_y, mode='lines', name='Fitted Curve', line=dict(color='rgba(0,100,80,1)', dash='solid', width=2)))

    plot_df = pd.DataFrame({'x': empirical_x, 'y': empirical_y, 'group': species_groups, 'species': species_names})
    for group_name, style in MARKER_STYLES.items():
        group_df = plot_df[plot_df['group'] == group_name]
        if not group_df.empty:
            fig.add_trace(go.Scatter(x=group_df['x'], y=group_df['y'], mode='markers', name=group_name, marker=dict(symbol=style['symbol'], color=style['color'], size=10, line=dict(width=1, color='DarkSlateGrey')), hovertext=[f"Group: {g}<br>Species: {sp}<br>Conc: {x:.3g} {unit}" for g, sp, x in zip(group_df['group'], group_df['species'], group_df['x'])], hoverinfo='text'))
    
    p_value_percent = plot_data['p_value'] * 100
    if hcp is not None and hcp > 0: fig.add_trace(go.Scatter(x=[hcp], y=[p_value_percent], mode='markers', marker=dict(color='gold', size=14, symbol='star', line=dict(color='black', width=1)), name=f'HC{p_value_percent:.0f}'))
    
    fig.update_layout(
        title=title, 
        xaxis_title=f'Concentration ({unit})', 
        yaxis_title='Percent of Species Affected (%)', 
        xaxis=dict(type="log", showgrid=False),
        yaxis=dict(range=[0, 100], showgrid=False), 
        legend_title='Legend', 
        template='plotly_white'
    )
    return fig

def render_diagnostics_table(results_df, hcp_percentile):
    st.markdown("This table shows the goodness-of-fit for the model(s) used in the analysis.")
    display_cols = ['name', 'aicc', 'weight', 'hcp', 'ks_pvalue', 'ad_statistic']
    display_df = results_df[display_cols].copy()
    display_df.rename(columns={'name': 'Distribution', 'aicc': 'AICc', 'weight': 'Weight', 'hcp': f'HC{hcp_percentile:.0f}', 'ks_pvalue': 'KS p-value', 'ad_statistic': 'AD Statistic'}, inplace=True)
    
    if MATPLOTLIB_INSTALLED:
        st.dataframe(display_df.style.format(precision=3).background_gradient(cmap='Greens', subset=['Weight']).highlight_min(subset=['AICc', 'AD Statistic'], color='#F9EBEA').highlight_max(subset=['KS p-value'], color='#E8F8F5'), use_container_width=True)
    else:
        st.warning("Matplotlib not installed. Displaying basic table. To enable color styling, run: `pip install matplotlib`")
        st.dataframe(display_df, use_container_width=True)
    return display_df