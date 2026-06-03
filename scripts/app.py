import os
import streamlit as st

st.set_page_config(
    page_title="FungiFlow",
    page_icon="🍄",
    layout="wide",
    initial_sidebar_state="expanded",
)

DATA_ROOT = "/data"


def load_css():
    css_path = os.path.join(DATA_ROOT, "css", "style.css")
    if os.path.exists(css_path):
        with open(css_path, encoding="utf-8") as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)


load_css()

from blast_app import show_blast_page
from iqtree_app import show_iqtree_page
from mafft_app import show_mafft_page

PAGES = {
    "Identification (BLAST)": "UNITE taxonomy & BLAST QC",
    "Consensuses (MAFFT)": "Alignment viewer",
    "Build IQ-Tree": "Phylogenetic trees",
}

with st.sidebar:
    st.markdown(
        """
        <div class="ff-brand">
            <p class="ff-brand-title">FungiFlow</p>
            <p class="ff-brand-sub">Nanopore ITS pipeline</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    page = st.radio("Module", list(PAGES.keys()), label_visibility="collapsed")
    st.markdown("---")
    st.caption(PAGES[page])
    st.markdown("---")

    raw_dir = os.path.join(DATA_ROOT, "raw_data")
    blast_dir = os.path.join(DATA_ROOT, "blast_results")
    n_raw = (
        len([f for f in os.listdir(raw_dir) if f.endswith(".fastq.gz")])
        if os.path.isdir(raw_dir)
        else 0
    )
    n_blast = (
        len([f for f in os.listdir(blast_dir) if f.endswith("_blast_summary.tsv")])
        if os.path.isdir(blast_dir)
        else 0
    )

    # Native metrics — avoids invisible **bold** in captions
    st.metric(label="Barcodes (raw_data)", value=n_raw)
    st.metric(label="BLAST summaries", value=n_blast)

if page == "Identification (BLAST)":
    show_blast_page()
elif page == "Consensuses (MAFFT)":
    show_mafft_page()
elif page == "Build IQ-Tree":
    show_iqtree_page()
