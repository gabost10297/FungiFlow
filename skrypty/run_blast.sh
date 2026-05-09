#!/bin/bash

# Zmienne
BLAST_DB="/data/baza_ITS/unite_blast_db"
THREADS=10

echo "ROZPOCZYNAM LOKALNE WYSZUKIWANIE BLAST"
echo "====================================================="

# Główny folder na wyniki
mkdir -p /data/wyniki_blast

# Pętla do każdego podfolderu z konsensusami
for sample_dir in /data/wyniki_spoa_konsensus/*; do
    if [ -d "$sample_dir" ]; then
        
        # Wyciągamy nazwę próbki
        SAMPLE=$(basename "$sample_dir")
        echo ">>> BLASTuję próbkę: ${SAMPLE} <<<"
        
        # Tworzymy plik wynikowy TSV dla danej próbki
        OUT_FILE="/data/wyniki_blast/${SAMPLE}_blast_summary.tsv"
        
        # Wpisujemy nagłówki do tabeli
        echo -e "Nazwa_Klastra\tID_Referencji\tZgodnosc_Procentowa\tDlugosc_Dopasowania\tPokrycie_Zapytania(%)\tE-value\tNazwa_Gatunku" > "$OUT_FILE"
        
        # Dla każdego pliku fasta w folderze tej próbki uruchamiamy BLAST
        for fasta in "$sample_dir"/*.fasta; do
            if [ -f "$fasta" ]; then
                
                # Uruchamiamy lokalnego blasta
                # -max_target_seqs 1 -> Zwraca tylko najlepsze dopasowanie
                # -outfmt "6 ..." -> Format tabelaryczny z wybranymi kolumnami
                blastn -query "$fasta" -db "$BLAST_DB" \
                       -outfmt "6 qseqid sseqid pident length qcovs evalue stitle" \
                       -max_target_seqs 5 -num_threads $THREADS >> "$OUT_FILE"
            fi
        done
        
        echo "Gotowe: ${SAMPLE}! Zapisano w ${OUT_FILE}"
    fi
done

echo "====================================================="
echo "WSZYSTKIE PRÓBKI ZIDENTYFIKOWANE!"