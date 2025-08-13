# Species Sensitivity Distribution (SSD) Generator

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://sstac-ssd-web-app.streamlit.app/)

A user-friendly web application for creating Species Sensitivity Distributions (SSDs), a key tool in ecological risk assessment.

![Application Screenshot](assets/screenshot.png)
*(Note: To add a screenshot, place an image named `screenshot.png` in a new `assets` folder.)*

## Overview

This application provides an intuitive interface for environmental scientists and toxicologists to perform robust SSD analysis. It can generate distributions, derive Hazard Concentrations (HCp), and visualize the results interactively. The app supports direct connections to a Supabase toxicology database or analysis of user-provided CSV files.

### Key Features

**Data Handling & Filtering**
* **Dynamic Data Sources:** Connect live to a Supabase database or upload a local CSV file.
* **Chemical & Media Filtering:** Easily select chemicals and filter data for Freshwater or Marine environments.
* **Intelligent Data Aggregation:** Handle multiple toxicity values for a single species by calculating the **Geometric Mean** (standard) or using the **Most Sensitive** value.
* **Robust Endpoint Cleaning:** Automatically identifies and excludes non-numeric or invalid endpoint data (e.g., 'NR', 'Not Reported') to ensure analysis quality.

**Powerful Statistical Analysis**
* **Multiple Distributions:** Utilizes four standard distributions: Log-Normal, Log-Logistic, Weibull, and Gamma.
* **Advanced Model Averaging:** Fits all four distributions and uses Akaike Information Criterion (AICc) to calculate model weights, producing a more robust, weighted-average result.
* **Single Distribution Mode:** Option to analyze data using a single, user-selected distribution.
* **Bootstrap Confidence Intervals:** Generates 95% confidence intervals for both the HCp value and the fitted distribution curve using a non-parametric bootstrap procedure.

**Interactive Outputs & UX**
* **Interactive Plotly Graph:** Visualizes the SSD curve, confidence bands, and individual data points styled by taxonomic group.
* **Customizable Protection Level (HCp):** The percentile (`p` in HCp) is fully customizable, allowing for the calculation of any hazard concentration (e.g., HC5, HC10, HC20).
* **Downloadable Results:** Export model diagnostic tables (GOF statistics) and the final aggregated data used for plotting.
* **Live Progress Bar:** A real-time status box provides feedback during computationally intensive bootstrap iterations.

## Core Scientific Methodology

The core analysis is powered by the `scipy` and `numpy` libraries, following established best practices.

1.  **Fitting Strategy:** To ensure proper statistical treatment, the application uses two fitting strategies:
    * **Log-Transformed Data:** The `Log-Normal` and `Log-Logistic` distributions are fitted to the natural logarithm of the toxicity data.
    * **Original Data:** The `Weibull` and `Gamma` distributions are fitted to the toxicity data in its original scale.

2.  **Model Parameterization:** All distributions are fitted using `scipy.stats.fit`, which estimates three parameters by default: shape, location (`loc`), and scale. The Akaike Information Criterion calculation correctly accounts for these parameters to ensure proper model penalization.

3.  **Model Averaging (AICc):** When model averaging is selected, the weight of each distribution is determined by its AICc score. This method rewards models for goodness-of-fit while penalizing them for complexity (number of parameters), preventing overfitting and providing a more reliable, averaged estimate.

## Project Structure

* `ssd_app.py`: The main Streamlit application file containing the UI and workflow logic.
* `ssd_core.py`: Handles the core statistical calculations, including distribution fitting and bootstrapping.
* `database.py`: Manages the connection to the Supabase database and data fetching.
* `ui_components.py`: Contains functions for generating the Plotly graph and other UI elements.
* `utils.py`: Includes helper functions, such as taxonomic group mapping and AICc calculation.

## Getting Started

### Prerequisites
* Python 3.9+
* An environment manager like `venv` or `conda`.

### Installation and Setup

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/jasennelson/ssd_web_app.git](https://github.com/jasennelson/ssd_web_app.git)
    cd ssd_web_app
    ```
2.  **Create and activate a virtual environment:**
    ```bash
    # For Windows
    python -m venv venv
    .\venv\Scripts\activate

    # For macOS/Linux
    python3 -m venv venv
    source venv/bin/activate
    ```
3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```
4.  **Set up Supabase credentials (Optional, for database mode):**
    * Create a file at `.streamlit/secrets.toml`.
    * Add your Supabase URL and Publishable key:
        ```toml
        [connections.supabase]
        url = "YOUR_SUPABASE_URL"
        key = "YOUR_SUPABASE_PUBLISHABLE_KEY"
        ```

### Running the Application

1.  **Run the Streamlit app from your terminal:**
    ```bash
    streamlit run ssd_app.py
    ```
2.  Open the local URL provided by Streamlit in your web browser.

## Contributing

Contributions, bug reports, and feature requests are welcome! Please feel free to open an issue or submit a pull request on the GitHub repository.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.