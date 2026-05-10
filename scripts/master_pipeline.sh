#!/bin/bash

SAMPLE="barcode25"
THREADS=10
KRAKEN_DB="/data/database"

echo "Starting TEST pipeline for sample: ${SAMPLE}"
echo "-------------------------------------------------"

mkdir -p /data/kraken_results
mkdir -p /data/consensus_results
mkdir -p /data/intermediate_data
mkdir -p /data/intermediate_data/clusters4_${SAMPLE}

# 1. Adapter trimming
echo "[1/6] Adapter trimming (Porechop)..."
porechop_abi --ab_initio -i /data/raw_data/${SAMPLE}.fastq.gz -o /data/intermediate_data/${SAMPLE}_porechop.fastq.gz -t $THREADS > /dev/null 2>&1

# 2. Conversion to FASTA + HYGIENE (Removing reads below 300bp)
echo "[2/6] Converting to FASTA and discarding short reads..."
seqtk seq -a -L 300 /data/intermediate_data/${SAMPLE}_porechop.fastq.gz | awk '/^>/{print ">"substr($0,2,18);next} {print}' > /data/intermediate_data/${SAMPLE}_short18.fasta

# 3. Clustering
echo "[3/6] Read clustering (CD-HIT) at 98%..."
cd-hit-est -i /data/intermediate_data/${SAMPLE}_short18.fasta -o /data/intermediate_data/${SAMPLE}_cd_hit3 -c 0.98 -n 10 -aS 0.8 -aL 0.8 -T $THREADS -M 0 > /dev/null 2>&1

echo "Unpacking clusters..."
make_multi_seq.pl /data/intermediate_data/${SAMPLE}_short18.fasta /data/intermediate_data/${SAMPLE}_cd_hit3.clstr /data/intermediate_data/clusters4_${SAMPLE}/ > /dev/null 2>&1

for plik in /data/intermediate_data/clusters4_${SAMPLE}/*; do
    if [ -f "$plik" ] && [[ "$plik" != *.fasta ]]; then
        mv "$plik" "${plik}.fasta"
    fi
done

# 4. Consensus generation
echo "[4/6] Generating final sequences (SPOA)..."
for f in /data/intermediate_data/clusters4_${SAMPLE}/*.fasta; do
    if [ -f "$f" ]; then
        base=$(basename "$f" .fasta)
        
        nreads=$(grep -c "^>" "$f" | tr -d '\r\n[:space:]')
        
        # If the variable is not empty and has 20 or more reads - consensus
        if [[ -n "$nreads" ]] && (( nreads >= 20 )); then
            spoa "$f" > /data/intermediate_data/temp.fasta
            sed -i "1s|^>.*|>cob_${base}_consensus|" /data/intermediate_data/temp.fasta
            mv /data/intermediate_data/temp.fasta "/data/consensus_results/${base}.fasta"
        fi
    fi
done

# 5. KRAKEN 
echo "[5/6] Taxonomic analysis (Kraken2)..."
kraken2 --db $KRAKEN_DB \
        --threads $THREADS \
        --report /data/kraken_results/${SAMPLE}_ITS_report.txt \
        --report-minimizer-data \
        --minimum-hit-groups 3 \
        /data/intermediate_data/${SAMPLE}_porechop.fastq.gz > /data/kraken_results/${SAMPLE}_ITS_kraken2.txt

# 6. CLEANUP 
echo "[6/6] Cleanup..."
# rm -rf /data/intermediate_data/clusters4_${SAMPLE}
# rm -f /data/intermediate_data/${SAMPLE}_*
# rm -f /data/intermediate_data/temp.fasta
rm -rf /data/tmp

echo "====================================================="
echo "Completed! Results in /kraken_results and /consensus_results"