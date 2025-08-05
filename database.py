# database.py
import streamlit as st
from st_supabase_connection import SupabaseConnection
import pandas as pd
import math

@st.cache_resource(show_spinner="Connecting to database...")
def initialize_supabase():
    """Initializes and returns the Supabase connection using Streamlit's secrets."""
    try:
        url = st.secrets.connections.supabase.url
        key = st.secrets.connections.supabase.key
        return st.connection("supabase", type=SupabaseConnection, url=url, key=key)
    except Exception as e:
        st.error(f"Failed to connect to Supabase. Ensure URL/Key are correct in secrets.toml. Error: {e}")
        return None

@st.cache_data(ttl=60)
def search_chemicals_in_db(_supabase_client, search_term: str):
    """Performs a two-stage search for relevance and returns a combined list."""
    if not _supabase_client or not search_term: return []
    try:
        prefix_response = _supabase_client.table("toxicology_data").select("chemical_name").ilike("chemical_name", f"{search_term}%").limit(50).execute()
        prefix_results = [item['chemical_name'] for item in prefix_response.data if item['chemical_name']]

        substring_response = _supabase_client.table("toxicology_data").select("chemical_name").ilike("chemical_name", f"%{search_term}%").limit(50).execute()
        substring_results = [item['chemical_name'] for item in substring_response.data if item['chemical_name']]

        combined_list = list(dict.fromkeys(prefix_results + substring_results))
        return sorted(combined_list)
            
    except Exception as e:
        st.error(f"Database search failed: {e}")
    return []

# --- DEFINITIVE FIX: HIGH-PERFORMANCE PAGINATION WITHOUT COUNTING ---
def fetch_data_for_chemicals(_supabase_client, chemical_names: list):
    """Fetches data for chemicals using robust pagination that does not require a pre-count."""
    if not chemical_names: return pd.DataFrame()

    required_columns = ['chemical_name', 'species_scientific_name', 'conc1_mean', 'species_group', 'media_type', 'endpoint', 'publication_year', 'author', 'title']
    all_data = []
    page = 0
    page_size = 1000  # Supabase's default and max limit per request

    # Create a persistent spinner
    spinner = st.spinner("Fetching data...")
    
    with spinner:
        while True:
            start_row = page * page_size
            end_row = start_row + page_size - 1
            
            try:
                # Update the spinner text for each page request
                spinner.text = f"Fetching data page {page + 1}..."
                
                page_response = _supabase_client.table("toxicology_data") \
                                                .select(','.join(required_columns)) \
                                                .in_("chemical_name", chemical_names) \
                                                .range(start_row, end_row) \
                                                .execute()
                
                if page_response.data:
                    all_data.extend(page_response.data)
                    # If we receive fewer rows than we asked for, it must be the last page.
                    if len(page_response.data) < page_size:
                        break
                else:
                    # If we receive no data, we are done.
                    break
                
                page += 1

            except Exception as e:
                st.error(f"Database fetch failed on page {page + 1}: {e}")
                return pd.DataFrame()

    return pd.DataFrame(all_data) if all_data else pd.DataFrame()