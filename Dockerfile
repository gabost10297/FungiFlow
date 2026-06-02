# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:1.5.10

USER root

RUN micromamba config append channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config set channel_priority strict

RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba install -n base -y \
      python=3.10 \
      porechop_abi cd-hit spoa seqtk gawk \
      kraken2 blast \
      mafft trimal iqtree \
      r-base r-ggplot2 bioconductor-ggtree bioconductor-ggtreeextra \
    && micromamba clean -a -y

RUN micromamba run -n base pip install --no-cache-dir \
      "streamlit>=1.28,<2" \
      "pandas>=2.0,<3" \
      "plotly>=5.0,<6"

WORKDIR /data
CMD ["bash"]