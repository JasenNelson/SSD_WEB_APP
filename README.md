# Species Sensitivity Distribution (SSD) Generator v3.0

This Streamlit web application generates Species Sensitivity Distribution (SSD) plots and calculates Hazard Concentrations (HCp) for aquatic life protection. Version 3.0 is a major overhaul, implementing more robust statistical techniques aligned with established practices, such as those found in the `ssdtools` R package.

The application allows users to analyze ecotoxicity data from either an uploaded CSV file or a built-in database powered by Supabase. It processes the data, fits multiple statistical distributions, and uses model averaging to derive a reliable water quality guideline.

## Key Features in Version 3.0

This version introduces significant methodological enhancements for more accurate and defensible results.

* [cite_start]**Model Averaging** [cite: 3][cite_start]: Instead of relying on a single user-selected distribution, the application now automatically fits four different distributions (Log-Normal, Log-Logistic, Gamma, and Weibull)[cite: 3]. [cite_start]It then calculates a model-averaged Hazard Concentration (HC5) based on Akaike Information Criterion (AICc) weights, providing a more robust estimate that accounts for model uncertainty[cite: 3].
* [cite_start]**Goodness-of-Fit (GoF) Statistics** [cite: 3][cite_start]: The application now performs Anderson-Darling goodness-of-fit tests for each distribution[cite: 3]. [cite_start]The results, including p-values and model weights, are presented in a clear table, allowing users to transparently assess how well each model fits the data[cite: 3].
* [cite_start]**Robust Guideline Calculation** [cite: 3][cite_start]: The final recommended guideline is now derived from the model-averaged 95% confidence interval, a more statistically sound approach[cite: 3]. [cite_start]The original "Protection Clause" (which ensures the guideline is no higher than the most sensitive species) is retained as a final conservative check on the output[cite: 3].
* [cite_start]**Corrected Data Filtering** [cite: 3][cite_start]: The data filtering logic for water types (Freshwater vs. Marine) has been fixed to ensure accurate subsetting of data based on user selection[cite: 3].
* **Interactive Data Visualization**: Generates high-quality, interactive SSD plots using Plotly. The plot displays the individual data points (color-coded by taxonomic group), the model-averaged fit line, and the 95% confidence interval.
* **Flexible Data Handling**: Offers options to aggregate multiple data points per species using either the geometric mean or the most sensitive (minimum) value.

## How to Use

All configuration options are located in the sidebar.

1.  **Choose Data Source**: Select "Upload File" to use your own CSV or "Database Search" to query the built-in library.
2.  **Select Chemical(s)**: Choose one or more chemicals to include in the analysis.
3.  **Set Guideline & Filter Options**:
    * **Water Type**: Filter data for Freshwater (FW), Marine (MW), or both.
4.  **Configure SSD Parameters**:
    * **Handle Multiple Values per Species**: Choose the aggregation method (Geometric Mean or Most Sensitive).
    * **Minimum number of species**: Set the minimum data points required (defaults to 8).
    * **Required Taxonomic Groups**: Select the taxonomic groups that must be present in the dataset.
5.  **Set Protection & Safety**:
    * **Hazard Concentration (HCp) Percentile**: Set the desired percentile (e.g., 5.0 for HC5).
6.  **Generate SSD**: Click the button to run the analysis.

## Methodology

### Data Aggregation and Filtering

The application first filters the dataset based on the selected chemical, water type, and exposure term. For long-term guidelines, it selects the most sensitive endpoint for each species based on a predefined preference ranking. It then aggregates multiple toxicity values for a single species using either the geometric mean or the minimum value, as specified by the user.

### Model Fitting and Averaging

The core of the analysis follows these steps:

1.  [cite_start]**Multi-Distribution Fitting**: The processed data is fit to four separate distributions: Log-Normal, Log-Logistic, Gamma, and Weibull[cite: 3].
2.  **Goodness-of-Fit Testing**: An Anderson-Darling test is performed for each fitted distribution to assess how well it describes the data.
3.  **AICc Weight Calculation**: The Akaike Information Criterion, corrected for small sample sizes (AICc), is calculated for each model. These values are then used to determine a weight for each model, where models that fit the data better receive higher weights.
4.  **Model-Averaged HCp**: The final HCp and its 95% confidence interval are calculated as a weighted average across all four models, using the AICc weights. This approach provides a more robust estimate than relying on a single distribution.

The final recommended guideline is the lower 95% confidence limit of the model-averaged HCp, ensuring a protective estimate.

## Requirements

The application is built on Python 3.9+ and requires the following key libraries:

* `streamlit`
* `pandas`
* `numpy`
* `scipy`
* `plotly`
* `statsmodels` (new in v3.0)
* `supabase`