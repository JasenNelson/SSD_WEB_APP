# Species Sensitivity Distribution (SSD) Generator

This web application, built with Streamlit, generates Species Sensitivity Distributions (SSDs) from toxicology data. It is designed to provide a user-friendly interface for environmental scientists and researchers to perform SSD analysis, visualize the results, and derive Hazard Concentrations (HCp) for ecological risk assessment.

The application can connect directly to a toxicology database hosted on Supabase or use a user-provided CSV file.

## How it Works

The application performs SSD analysis by fitting one or more statistical distributions to the provided toxicity data. The core analysis is powered by the `scipy.stats` library.

* **Distributions:** It supports four common distributions used in ecotoxicology: Log-Normal, Log-Logistic, Weibull, and Gamma.
* **Model Fitting:** The application fits the selected distribution(s) to the species toxicity data to generate a cumulative distribution function (CDF).
* **Model Averaging:** In "Model Averaging" mode, the application fits all four distributions and calculates weights based on the Akaike Information Criterion corrected for small sample sizes (AICc). The final HCp and CDF curve are a weighted average of all models, providing a more robust estimate.
* **Confidence Intervals:** 95% confidence intervals for the HCp value and the fitted curve are generated using a non-parametric bootstrapping method.

## Key Features

* **Dynamic Data Sourcing:**
    * Connect directly to a Supabase toxicology database with a live chemical search.
    * Upload local data via CSV file for analysis.
* **Advanced Filtering:**
    * Filter data by chemical and water type (Freshwater or Marine).
    * Robustly excludes data with "Not Reported" endpoints (e.g., 'NR', 'NR-ZERO') to ensure analysis quality.
* **Customizable Analysis:**
    * Handle multiple toxicity values per species by taking the **Geometric Mean** or the **Most Sensitive (Minimum)** value.
    * Choose between a robust **Model Averaging** approach or analysis with a **Single Distribution**.
* **Interactive Visualization & Output:**
    * Displays an interactive SSD plot using Plotly, with data points categorized and styled by taxonomic group.
    * Calculates and displays the HCp estimate and its 95% confidence interval.
    * Provides downloadable tables for model diagnostics and the final aggregated data used in the plot.
* **Enhanced User Experience:**
    * A **"Select All"** checkbox to easily analyze all available chemicals.
    * A detailed status box with a **real-time progress bar** provides feedback during the bootstrap analysis.

## Project Structure

The application is organized into several key Python files:

* `ssd_app.py`: The main Streamlit application file containing the UI and workflow logic.
* `core_analysis.py`: Handles the core statistical calculations for the SSD, including distribution fitting and bootstrapping.
* `database.py`: Manages the connection to the Supabase database and data fetching.
* `ui_components.py`: Contains functions for generating the Plotly graph and other UI elements.
* `utils.py`: Includes helper functions, such as the taxonomic group mapping.

## Installation and Setup

1.  **Clone the repository:**
    ```bash
    git clone <your-repository-url>
    cd <your-repository-folder>
    ```
2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```
3.  **Set up Supabase credentials (Required for Database Search):**
    * Create a file at `.streamlit/secrets.toml`.
    * Add your Supabase URL and Key to this file:
        ```toml
        [connections.supabase]
        url = "YOUR_SUPABASE_URL"
        key = "YOUR_SUPABASE_ANON_KEY"
        ```

## Usage

1.  **Run the Streamlit app:**
    ```bash
    streamlit run ssd_app.py
    ```
2.  In the sidebar, choose your **Data Source**.
3.  **Select Chemicals** for analysis. You can use the search bar (for database mode) and the "Select All" checkbox.
4.  Configure **Filtering** (Water Type), **SSD Parameters** (Aggregation, Analysis Mode), and **Protection** level (HCp Percentile, Bootstrap Iterations).
5.  Click the **"Generate SSD"** button.
6.  Monitor the progress in the status box that appears.
7.  Once complete, review your results in the main panel, which is organized into four tabs:
    * **Summary & Plot:** The main results and interactive SSD graph.
    * **Model Diagnostics:** Goodness-of-fit statistics for the models.
    * **Final Data:** The aggregated data used for plotting.
    * **Processing Log:** Any warnings generated during the analysis.

## Data Requirements

For the app to function correctly, the input data (either from the database or CSV) must contain the following columns:
* `chemical_name` (text)
* `species_scientific_name` (text)
* `conc1_mean` (numeric) - The toxicity value.
* `species_group` (text) - e.g., 'Fish', 'Plant', 'Invertebrate'.
* `endpoint` (text)
* `media_type` (text) - Should contain 'FW' for freshwater or 'MW' for marine.