import streamlit as st
import os

# Dahboard page configuration
st.set_page_config(page_title="FungiFlow Dashboard", layout="wide")

from blast_app import show_blast_page
from mafft_app import show_mafft_page

st.sidebar.title("Navigation")
st.sidebar.markdown("---")

# Module choice
page = st.sidebar.radio(
    "Choose analysis module: ", 
    ["Identification (BLAST)", "Consensuses (MAFFT)"]
)
st.sidebar.markdown("---")

if page == "Identification (BLAST)":
    show_blast_page()
elif page == "Consensuses (MAFFT)":
    show_mafft_page()

def load_css():
    css_path = "/data/css/style.css"
    if os.path.exists(css_path):
        with open(css_path) as f:
            st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

load_css()