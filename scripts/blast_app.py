import streamlit as st
import pandas as pd
import plotly.express as px
import glob
import os
import subprocess

# --- FUNKCJE POMOCNICZE ---
@st.cache_data
def load_and_clean_data(path):
    df = pd.read_csv(path, sep='\t')
    def parse_taxonomy(tax_string):
        try:
            parts = str(tax_string).split('|')
            tax_raw = parts[-1] if len(parts) > 1 else str(tax_string)
            ranks = tax_raw.split(';')
            tax_dict = {}
            for r in ranks:
                if '__' in r:
                    lvl, val = r.split('__', 1)
                    tax_dict[lvl] = val
            return pd.Series(tax_dict)
        except:
            return pd.Series({})

    if 'Species_Name' in df.columns:
        parsed_tax = df['Species_Name'].apply(parse_taxonomy)
        parsed_tax = parsed_tax.rename(columns={
            'k': 'Kingdom', 'p': 'Phylum', 'c': 'Class', 
            'o': 'Order', 'f': 'Family', 'g': 'Genus', 's': 'Species'
        })
        df_clean = pd.concat([df, parsed_tax], axis=1)
    else:
        df_clean = df
    return df_clean.drop_duplicates(subset=['Cluster_Name'], keep='first')

def show_blast_page():
    st.title("BLAST Results Explorer")
    st.info("Taxonomic analysis and quality metrics for individual samples (barcodes).")

    # Ścieżki
    blast_log_path = "/data/blast_results/blast_run.log"
    blast_done_flag = "/data/blast_results/blast_done.flag"
    
    with st.expander("BLAST Execution Settings"):
        if os.path.exists(blast_log_path) and not os.path.exists(blast_done_flag):
            st.warning("BLAST is currently running in the background...")
            if st.button("Refresh Status"): st.rerun()
            try:
                with open(blast_log_path, "r") as f:
                    st.code("".join(f.readlines()[-10:]))
            except: pass
        else:
            if st.button("Run BLAST Analysis", type="primary", key="blast_run_button_unique"):
                with open(blast_log_path, "w") as log_file:
                    # Wywołanie skryptu z Twojej instrukcji:
                    subprocess.run(["bash", "/data/scripts/run_blast.sh"], stdout=log_file, stderr=subprocess.STDOUT)
                with open(blast_done_flag, "w") as f: f.write("done")
                st.rerun()

    files = glob.glob("/data/blast_results/*_blast_summary.tsv")
    if not files:
        st.info("No results found. Run BLAST Analysis above.")
        return

    selected_file = st.selectbox("Select sample to analyze:", [os.path.basename(f) for f in files])
    df = load_and_clean_data(f"/data/blast_results/{selected_file}")

    col1, col2 = st.columns(2)
    with col1:
        if 'Genus' in df.columns:
            st.subheader("Genus Distribution")
            fig = px.pie(df['Genus'].value_counts().reset_index(), values='count', names='Genus', hole=0.3)
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Quality Metrics")
        if 'Percent_Identity' in df.columns:
            fig = px.scatter(df, x='Percent_Identity', y='Alignment_Length', color='Genus')
            fig.update_xaxes(autorange="reversed")
            st.plotly_chart(fig, use_container_width=True)

    st.subheader("Top 1 Matches")
    st.dataframe(df.reset_index(drop=True), use_container_width=True)

show_blast_page()