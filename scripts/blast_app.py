import streamlit as st
import pandas as pd
import plotly.express as px
import glob
import os

# Function to clean data (UNITE taxonomy)
@st.cache_data
def load_and_clean_data(path):
    df = pd.read_csv(path, sep='\t')
    
    # Function extracting specific ranks from UNITE name
    def parse_taxonomy(tax_string):
        try:
            # Split by "|" and take the last element
            parts = str(tax_string).split('|')
            tax_raw = parts[-1] if len(parts) > 1 else str(tax_string)
            
            # Split by ";" and map ranks
            ranks = tax_raw.split(';')
            tax_dict = {}
            for r in ranks:
                if '__' in r:
                    lvl, val = r.split('__', 1)
                    tax_dict[lvl] = val
            return pd.Series(tax_dict)
        except:
            return pd.Series({})

    # Appending new taxonomy columns to the table
    if 'Species_Name' in df.columns:
        parsed_tax = df['Species_Name'].apply(parse_taxonomy)
        parsed_tax = parsed_tax.rename(columns={
            'k': 'Kingdom', 'p': 'Phylum', 'c': 'Class', 
            'o': 'Order', 'f': 'Family', 'g': 'Genus', 's': 'Species'
        })
        df_clean = pd.concat([df, parsed_tax], axis=1)
    else:
        df_clean = df
    
    # Keep only the first result for each cluster
    df_top1 = df_clean.drop_duplicates(subset=['Cluster_Name'], keep='first')
    return df_top1


def show_blast_page():
    st.title("BLAST Results Explorer")
    st.info("Taxonomic analysis and quality metrics for individual samples (barcodes).")

    # Find result files
    files = glob.glob("/data/blast_results/*_blast_summary.tsv")
    if not files:
        st.warning("No TSV files found in the /data/blast_results/ folder.")
        return # Return instead of st.stop() in a module

    st.sidebar.subheader("BLAST Options")
    # Select sample from the dropdown list in the sidebar
    selected_file = st.sidebar.selectbox("Select sample to analyze:", [os.path.basename(f) for f in files])
    file_path = f"/data/blast_results/{selected_file}"

    df = load_and_clean_data(file_path)

    # DATA VISUALIZATION
    st.markdown("---")
    col1, col2 = st.columns(2)

    with col1:
        with st.container(border=True):
            st.subheader("Genus Distribution")
            if 'Genus' in df.columns:
                genus_counts = df['Genus'].value_counts().reset_index()
                genus_counts.columns = ['Genus', 'Cluster count']
                
                fig_pie = px.pie(genus_counts, values='Cluster count', names='Genus', hole=0.3)
                
                fig_pie.update_layout(
                    legend=dict(orientation="h", yanchor="top", y=-0.1, xanchor="center", x=0.5),
                    margin=dict(t=20, b=20, l=20, r=20)
                )
                st.plotly_chart(fig_pie, use_container_width=True)

    with col2:
        with st.container(border=True):
            st.subheader("Quality Metrics")
            if 'Percent_Identity' in df.columns and 'Alignment_Length' in df.columns:
                fig_scatter = px.scatter(df, x='Percent_Identity', y='Alignment_Length', 
                                         color='Genus' if 'Genus' in df.columns else None, 
                                         hover_data=['Species', 'Cluster_Name'] if 'Species' in df.columns else None)
                # Reverse X axis
                fig_scatter.update_xaxes(autorange="reversed")
                fig_scatter.update_layout(
                    legend=dict(orientation="h", yanchor="top", y=-0.2, xanchor="center", x=0.5),
                    margin=dict(t=20, b=20, l=20, r=20)
                )
                st.plotly_chart(fig_scatter, use_container_width=True)

    st.markdown("---")
    with st.container(border=True):
        st.subheader("Top 1 Match")
        
        # TABLE COSMETICS
        df_display = df.copy()
        
        columns_to_drop = ['Reference_ID', 'Species_Name']
        df_display = df_display.drop(columns=[k for k in columns_to_drop if k in df_display.columns])
        
        df_display = df_display.reset_index(drop=True)
        df_display.index = df_display.index + 1

        if 'E-value' in df_display.columns:
            df_display['E-value'] = pd.to_numeric(df_display['E-value'], errors='coerce')
            df_display['E-value'] = df_display['E-value'].apply(lambda x: f"{x:.2e}" if x != 0 else "0.0")
        
        st.dataframe(df_display, use_container_width=True)