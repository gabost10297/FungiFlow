#!/bin/bash

# Variables
BLAST_DB="/data/database/unite_blast_db"
THREADS=10

echo "STARTING LOCAL BLAST SEARCH"
echo "====================================================="

# Main output directory
mkdir -p /data/blast_results

# Loop through each consensus folder
for sample_dir in /data/consensus_results/*; do
    if [ -d "$sample_dir" ]; then
        
        # Extract the sample name
        SAMPLE=$(basename "$sample_dir")
        echo ">>> BLASTing sample: ${SAMPLE} <<<"
        
        # Create TSV output file for the given sample
        OUT_FILE="/data/blast_results/${SAMPLE}_blast_summary.tsv"
        
        # Write headers to the table
        echo -e "Cluster_Name\tReference_ID\tPercent_Identity\tAlignment_Length\tQuery_Coverage(%)\tE-value\tSpecies_Name" > "$OUT_FILE"
        
        # Run BLAST for each fasta file in the sample folder
        for fasta in "$sample_dir"/*.fasta; do
            if [ -f "$fasta" ]; then
                
                # Run local BLAST
                # -max_target_seqs 1 -> Returns only the best match
                # -outfmt "6 ..." -> Tabular format with selected columns
                # 2> /dev/null -> Hides warnings
                blastn -query "$fasta" -db "$BLAST_DB" \
                       -outfmt "6 qseqid sseqid pident length qcovs evalue stitle" \
                       -max_target_seqs 1 -num_threads $THREADS >> "$OUT_FILE" 2> /dev/null
            fi
        done
        
        echo "Done: ${SAMPLE}! Saved to ${OUT_FILE}"
    fi
done

echo "====================================================="
echo "ALL SAMPLES IDENTIFIED!"