# Baza z Condą
FROM continuumio/miniconda3:latest

# Dodajemy oficjalne sklepy
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Instalujemy narzędzia + WYMUSZAMY starszą wersję Pythona
RUN conda install -y python=3.10 \
    porechop_abi \
    cd-hit \
    spoa \
    vsearch \
    kraken2 \
    gawk \
    seqtk \
    && conda clean -afy

# Ustawiamy katalog roboczy
WORKDIR /data

# Domyślne zachowanie
CMD ["bash"]

# Biblioteki python
RUN pip install streamlit pandas plotly