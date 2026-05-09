#!/bin/bash

THREADS=10
KRAKEN_DB="/data/baza_ITS"

echo "ROZPOCZYNAM PRZETWARZANIE PRÓBEK"
echo "====================================================="

# Pętla po plikach .fastq.gz w folderze dane_surowe
for plik_full in /data/dane_surowe/*.fastq.gz; do
    
    # Wyciągnij samą nazwę
    SAMPLE=$(basename "$plik_full" .fastq.gz)
    
    echo ">>> ROZPOCZYNAM PRACĘ NAD PRÓBKĄ: ${SAMPLE} <<<"
    
    # Osobne podfoldery dla konkretnej próbki
    mkdir -p /data/wyniki_kraken/${SAMPLE}
    mkdir -p /data/wyniki_spoa_konsensus/${SAMPLE}
    mkdir -p /data/dane_posrednie
    mkdir -p /data/dane_posrednie/clusters4_${SAMPLE}

    # 1. Wycinki adapterow
    echo "[1/6] Wycinka adapterow (Porechop)..."
    porechop_abi --ab_initio -i "$plik_full" -o "/data/dane_posrednie/${SAMPLE}_porechop.fastq.gz" -t $THREADS > /dev/null 2>&1

    # 2. Konwersja do fasta i usuwanie poniżej 300
    echo "[2/6] Konwersja do FASTA i odrzucanie krotkich odczytow..."
    seqtk seq -a -L 300 "/data/dane_posrednie/${SAMPLE}_porechop.fastq.gz" | awk '/^>/{print ">"substr($0,2,18);next} {print}' > "/data/dane_posrednie/${SAMPLE}_short18.fasta"

    # 3. Klastrowanie
    echo "[3/6] Klastrowanie odczytow (CD-HIT) na poziomie 98%..."
    cd-hit-est -i "/data/dane_posrednie/${SAMPLE}_short18.fasta" -o "/data/dane_posrednie/${SAMPLE}_cd_hit3" -c 0.98 -n 10 -aS 0.8 -aL 0.8 -T $THREADS -M 0 > /dev/null 2>&1

    echo "Rozpakowywanie klastrów..."
    make_multi_seq.pl "/data/dane_posrednie/${SAMPLE}_short18.fasta" "/data/dane_posrednie/${SAMPLE}_cd_hit3.clstr" "/data/dane_posrednie/clusters4_${SAMPLE}/" > /dev/null 2>&1

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
            
            # ZLICZANIE ODCZYTÓW 
            nreads=$(grep -c "^>" "$f" | tr -d '\r\n[:space:]')
            
            # Jeśli zmienna nie jest pusta i ma 20 lub więcej odczytów - konsensus
            if [[ -n "$nreads" ]] && (( nreads >= 20 )); then
                spoa "$f" > "/data/dane_posrednie/temp_${SAMPLE}.fasta"
                # Podpisuje nagłówek nazwą próbki
                sed -i "1s|^>.*|>cob_${SAMPLE}_${base}_consensus|" "/data/dane_posrednie/temp_${SAMPLE}.fasta"
                # Wrzuca do podfolderu
                mv "/data/dane_posrednie/temp_${SAMPLE}.fasta" "/data/wyniki_spoa_konsensus/${SAMPLE}/${base}.fasta"
            fi
        fi
    done

    # 5. KRAKEN 
    echo "[5/6] Analiza taksonomiczna (Kraken2)..."
    kraken2 --db $KRAKEN_DB \
            --threads $THREADS \
            --report "/data/wyniki_kraken/${SAMPLE}/${SAMPLE}_ITS_report.txt" \
            --report-minimizer-data \
            --minimum-hit-groups 3 \
            "/data/dane_posrednie/${SAMPLE}_porechop.fastq.gz" > "/data/wyniki_kraken/${SAMPLE}/${SAMPLE}_ITS_kraken2.txt"

    # 6. SPRZĄTANIE
    echo "[6/6] Sprzatanie plikow tymczasowych dla ${SAMPLE}..."
    rm -rf "/data/dane_posrednie/clusters4_${SAMPLE}"
    rm -f "/data/dane_posrednie/${SAMPLE}_"*
    rm -f "/data/dane_posrednie/temp_${SAMPLE}.fasta"
    rm -rf /data/tmp

    echo "Zakonczono próbkę: ${SAMPLE}!"
    echo "====================================================="
done

echo "WSZYSTKIE PRÓBKI PRZETWORZONE POMYŚLNIE!"