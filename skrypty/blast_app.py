import streamlit as st
import pandas as pd
import plotly.express as px
import glob
import os

# Konfiguracja strony
st.set_page_config(page_title="Grzybowy Eksplorator", layout="wide")
st.title("Eksplorator Wyników BLAST")

# Znalezienie plików z wynikami
files = glob.glob("/data/wyniki_blast/*_blast_summary.tsv")
if not files:
    st.warning("Nie znaleziono plików TSV w folderze /data/wyniki_blast/")
    st.stop()

# Wybór próbki z rozwijanej listy
selected_file = st.selectbox("Wybierz próbkę do analizy:", [os.path.basename(f) for f in files])
file_path = f"/data/wyniki_blast/{selected_file}"

# Funkcja do czyszczenia danych (taksonomie z UNITE)
@st.cache_data
def load_and_clean_data(path):
    df = pd.read_csv(path, sep='\t')
    
    # Funkcja wyciągająca poszczególne rangi z nazwy UNITE
    def parse_taxonomy(tax_string):
        try:
            # Dzieli po znaku "|" i bierze ostatni element
            parts = str(tax_string).split('|')
            tax_raw = parts[-1] if len(parts) > 1 else str(tax_string)
            
            # Dzieli po ";" i mapuje rangi
            ranks = tax_raw.split(';')
            tax_dict = {}
            for r in ranks:
                if '__' in r:
                    lvl, val = r.split('__', 1)
                    tax_dict[lvl] = val
            return pd.Series(tax_dict)
        except:
            return pd.Series({})

    # Doklejanie nowych kolumn do tabeli
    parsed_tax = df['Nazwa_Gatunku'].apply(parse_taxonomy)
    parsed_tax = parsed_tax.rename(columns={
        'k': 'Królestwo', 'p': 'Typ', 'c': 'Klasa', 
        'o': 'Rząd', 'f': 'Rodzina', 'g': 'Rodzaj', 's': 'Gatunek'
    })
    
    df_clean = pd.concat([df, parsed_tax], axis=1)
    
    # Pierwszy z 5 wynikow dla kazdego klastra
    df_top1 = df_clean.drop_duplicates(subset=['Nazwa_Klastra'], keep='first')
    return df_top1

df = load_and_clean_data(file_path)

# WIZUALIZACJA DANYCH
st.markdown("---")
col1, col2 = st.columns(2)

with col1:
   
    with st.container(border=True):
        st.subheader(f"Rozkład Rodzajów")
        if 'Rodzaj' in df.columns:
            rodzaje_counts = df['Rodzaj'].value_counts().reset_index()
            rodzaje_counts.columns = ['Rodzaj', 'Liczba klastrów']
            
            fig_pie = px.pie(rodzaje_counts, values='Liczba klastrów', names='Rodzaj', hole=0.3)
           
            fig_pie.update_layout(
                legend=dict(orientation="h", yanchor="top", y=-0.1, xanchor="center", x=0.5),
                margin=dict(t=20, b=20, l=20, r=20)
            )
            st.plotly_chart(fig_pie, use_container_width=True)

with col2:
    with st.container(border=True):
        st.subheader("Wskaźniki Jakości")
        fig_scatter = px.scatter(df, x='Zgodnosc_Procentowa', y='Dlugosc_Dopasowania', 
                                 color='Rodzaj', hover_data=['Gatunek', 'Nazwa_Klastra'])
        # Odwrocenie osi X
        fig_scatter.update_xaxes(autorange="reversed")
        fig_scatter.update_layout(
            legend=dict(orientation="h", yanchor="top", y=-0.2, xanchor="center", x=0.5),
            margin=dict(t=20, b=20, l=20, r=20)
        )
        st.plotly_chart(fig_scatter, use_container_width=True)

st.markdown("---")
with st.container(border=True):
    st.subheader("Top 1 dopasowanie")
    
    # KOSMETYKA TABELI
    df_display = df.copy()
    
    kolumny_do_usuniecia = ['ID_Referencji', 'Nazwa_Gatunku']
    df_display = df_display.drop(columns=[k for k in kolumny_do_usuniecia if k in df_display.columns])
    
    df_display = df_display.reset_index(drop=True)
    df_display.index = df_display.index + 1

    if 'E-value' in df_display.columns:
        df_display['E-value'] = pd.to_numeric(df_display['E-value'], errors='coerce')
        df_display['E-value'] = df_display['E-value'].apply(lambda x: f"{x:.2e}" if x != 0 else "0.0")
    
    st.dataframe(df_display, use_container_width=True)