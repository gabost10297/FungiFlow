#!/bin/bash
set -e

# PATH CONFIGURATION - Inside the container
DB_DIR="/data/database"

echo "=== STARTING UNITE DATABASE BUILD ==="

# 0. CLEANUP
echo "Cleaning up old temporary files..."
rm -rf "$DB_DIR/library"
rm -rf "$DB_DIR/taxonomy"
rm -f "$DB_DIR"/*.k2d
rm -f "$DB_DIR"/unite_kraken.fasta
rm -f "$DB_DIR"/seqid2taxid.map
find "$DB_DIR" -name "*.fasta" -delete

# 1. EXTRACTING ARCHIVE
echo "[1/4] Extracting UNITE archive..."
TGZ_FILE=$(find "$DB_DIR" -maxdepth 1 -name "*.tgz" -o -name "*.tar.gz" | head -n 1)

if [ -z "$TGZ_FILE" ]; then
    echo "ERROR: No .tgz or .tar.gz file found in $DB_DIR!"
    exit 1
fi

echo "Found archive: $(basename "$TGZ_FILE")"
tar -xzf "$TGZ_FILE" -C "$DB_DIR"

# Find the extracted FASTA
RAW_FASTA=$(find "$DB_DIR" -name "*.fasta" | head -n 1)
FILE_NAME=$(basename "$RAW_FASTA")

# 2. CONVERTING (Now running directly, no 'docker run' needed!)
echo "[2/4] Running converter..."
python /data/scripts/unite_to_kraken.py "/data/database/$FILE_NAME"

# 3. ADDING TO LIBRARY (Direct call)
echo "[3/4] Adding sequences to Kraken library..."
kraken2-build --add-to-library /data/database/unite_kraken.fasta \
              --db /data/database/ \
              --no-masking

# 4. BUILDING DATABASE (Direct call)
echo "[4/4] Building Kraken dictionary..."
kraken2-build --build --db /data/database/ --threads 10

echo ""
echo "=== PROCESS COMPLETED ==="