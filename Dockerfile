# Base image with Conda
FROM continuumio/miniconda3:latest

# Add official channels (bioconda and conda-forge)
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Install bioinformatics tools (DODANO trimal)
RUN conda install -y python=3.10 \
    porechop_abi \
    cd-hit \
    spoa \
    vsearch \
    kraken2 \
    gawk \
    seqtk \
    mafft \
    trimal \
    iqtree \
    r-base \
    r-ggplot2 \
    bioconductor-ggtree \
    bioconductor-ggtreeextra \
    && conda clean -afy

# Install Python libraries for the Streamlit dashboard
RUN pip install --no-cache-dir streamlit pandas plotly

WORKDIR /data

CMD ["bash"]