#!/bin/bash
set -e

# KONFIGURACJA ŚCIEŻEK
PROJECT_DIR="/mnt/c/Users/zosra/OneDrive/Desktop/Projekt_Inzynierski"
DB_DIR="$PROJECT_DIR/baza_ITS"

echo "=== ROZPOCZYNAM BUDOWĘ BAZY UNITE ==="

# 0. CZYSZCZENIE
echo "Czyszczenie starych plików tymczasowych..."
rm -rf "$DB_DIR/library"
rm -rf "$DB_DIR/taxonomy"
rm -f "$DB_DIR"/*.k2d
rm -f "$DB_DIR"/unite_kraken.fasta
rm -f "$DB_DIR"/seqid2taxid.map
# Usuwamy stare wypakowane pliki FASTA
find "$DB_DIR" -name "*.fasta" -delete

# 1. ROZPAKOWANIE ARCHIWUM
echo "[1/4] Rozpakowywanie archiwum UNITE..."
tar -xzf "$DB_DIR/sh_general_release_19.02.2025.tgz" -C "$DB_DIR"

# Szukamy ścieżki do wypakowanego pliku FASTA
SUROWE_FASTA=$(find "$DB_DIR" -name "*.fasta" | head -n 1)
NAZWA_PLIKU=$(basename "$SUROWE_FASTA")

# 2. TŁUMACZENIE
echo "[2/4] Uruchamianie translatora..."
docker run --rm -v "${PROJECT_DIR}:/data" grzyby_pro \
    python /data/skrypty/unite_to_kraken.py "/data/baza_ITS/$NAZWA_PLIKU"

# 3. DODAWANIE DO BIBLIOTEKI
echo "[3/4] Wgrywanie sekwencji do biblioteki Krakena..."
docker run --rm -v "${PROJECT_DIR}:/data" grzyby_pro \
    kraken2-build --add-to-library /data/baza_ITS/unite_kraken.fasta \
                  --db /data/baza_ITS/ \
                  --no-masking

# 4. KOMPILACJA BAZY
echo "[4/4] Kompilacja słownika..."
docker run --rm -v "${PROJECT_DIR}:/data" grzyby_pro \
    kraken2-build --build --db /data/baza_ITS/ --threads 10

echo ""
echo "=== PROCES ZAKOŃCZONY ==="
echo "Weryfikacja: Jeśli w logach powyżej (Krok 4) widnieje liczba przetworzonych sekwencji"
echo "większa niż 0 (np. Completed processing of 147735 sequences), baza jest GOTOWA."