#!/bin/bash

INPUT_DIR="/data/consensus_results"
TMP_FASTA="/data/tmp/all_consensuses.fasta"
OUTPUT_FASTA="/data/intermediate_data/mafft_alignment.fasta"

echo "[1/2] Gathering consensus sequences..."
cat ${INPUT_DIR}/*/*.fasta > ${TMP_FASTA} 2>/dev/null || cat ${INPUT_DIR}/*.fasta > ${TMP_FASTA}

if [ ! -s "${TMP_FASTA}" ]; then
    echo "Error: No FASTA files containing consensus sequences were found in ${INPUT_DIR}!"
    exit 1
fi

echo "[2/2] Lunching comparison (MAFFT)..."
mafft --auto ${TMP_FASTA} > ${OUTPUT_FASTA}

echo "The result has been saved in: ${OUTPUT_FASTA}"