#!/bin/bash

SAMPLE="barcode25"
THREADS=10
KRAKEN_DB="/data/baza_ITS"

echo "Rozpoczynam TESTOWY potok dla próbki: ${SAMPLE}"
echo "-------------------------------------------------"

mkdir -p /data/wyniki_kraken
mkdir -p /data/wyniki_spoa_konsensus
mkdir -p /data/dane_posrednie
mkdir -p /data/dane_posrednie/clusters4_${SAMPLE}

# 1. Wycinki adapterow
echo "[1/6] Wycinka adapterow (Porechop)..."
porechop_abi --ab_initio -i /data/dane_surowe/${SAMPLE}.fastq.gz -o /data/dane_posrednie/${SAMPLE}_porechop.fastq.gz -t $THREADS > /dev/null 2>&1

# 2. Konwersja do fasta + HIGIENA (Usuwanie śmieci poniżej 300
echo "[2/6] Konwersja do FASTA i odrzucanie krotkich odczytow..."
seqtk seq -a -L 300 /data/dane_posrednie/${SAMPLE}_porechop.fastq.gz | awk '/^>/{print ">"substr($0,2,18);next} {print}' > /data/dane_posrednie/${SAMPLE}_short18.fasta

# 3. Klastrowanie
echo "[3/6] Klastrowanie odczytow (CD-HIT) na poziomie 98%..."
cd-hit-est -i /data/dane_posrednie/${SAMPLE}_short18.fasta -o /data/dane_posrednie/${SAMPLE}_cd_hit3 -c 0.98 -n 10 -aS 0.8 -aL 0.8 -T $THREADS -M 0 > /dev/null 2>&1

echo "Rozpakowywanie klastrów..."
make_multi_seq.pl /data/dane_posrednie/${SAMPLE}_short18.fasta /data/dane_posrednie/${SAMPLE}_cd_hit3.clstr /data/dane_posrednie/clusters4_${SAMPLE}/ > /dev/null 2>&1

for plik in /data/dane_posrednie/clusters4_${SAMPLE}/*; do
    if [ -f "$plik" ] && [[ "$plik" != *.fasta ]]; then
        mv "$plik" "${plik}.fasta"
    fi
done

# 4. Generowanie konsensusu
echo "[4/6] Generowanie ostatecznych sekwencji (SPOA)..."
for f in /data/dane_posrednie/clusters4_${SAMPLE}/*.fasta; do
    if [ -f "$f" ]; then
        base=$(basename "$f" .fasta)
        
        nreads=$(grep -c "^>" "$f" | tr -d '\r\n[:space:]')
        
        # Jeśli zmienna nie jest pusta i ma 20 lub więcej odczytów - konsensus
        if [[ -n "$nreads" ]] && (( nreads >= 20 )); then
            spoa "$f" > /data/dane_posrednie/temp.fasta
            sed -i "1s|^>.*|>cob_${base}_consensus|" /data/dane_posrednie/temp.fasta
            mv /data/dane_posrednie/temp.fasta "/data/wyniki_spoa_konsensus/${base}.fasta"
        fi
    fi
done

# 5. KRAKEN 
echo "[5/6] Analiza taksonomiczna (Kraken2)..."
kraken2 --db $KRAKEN_DB \
        --threads $THREADS \
        --report /data/wyniki_kraken/${SAMPLE}_ITS_report.txt \
        --report-minimizer-data \
        --minimum-hit-groups 3 \
        /data/dane_posrednie/${SAMPLE}_porechop.fastq.gz > /data/wyniki_kraken/${SAMPLE}_ITS_kraken2.txt

# 6. SPRZĄTANIE 
echo "[6/6] Sprzatanie..."
# rm -rf /data/dane_posrednie/clusters4_${SAMPLE}
# rm -f /data/dane_posrednie/${SAMPLE}_*
# rm -f /data/dane_posrednie/temp.fasta
rm -rf /data/tmp

echo "====================================================="
echo "Zakonczono! Wyniki w /wyniki_kraken oraz /wyniki_spoa_konsensus"