#  Potok Analizy ITS Grzybów (ITS Nanopore)

Zautomatyzowany system do przetwarzania surowych odczytów z sekwencjonowania Nanopore, generowania wysokiej jakości sekwencji konsensusowych oraz ich precyzyjnej identyfikacji taksonomicznej przy użyciu referencyjnej bazy UNITE. Projekt zawiera interaktywny dashboard do wizualizacji i analizy wyników.

## Struktura Projektu

Zgodnie z obecną konfiguracją, projekt wykorzystuje poniższy układ katalogów:

```text
PROJEKT_INZYNIERSKI/
├── .gitignore                  # Plik wykluczający duże dane z repozytorium GitHub
├── Dockerfile                  # Definicja obrazu kontenera (biotools + Python/Streamlit)
├── ITS_list.txt                # Lista pomocnicza
├── baza_ITS/                   # Bazy UNITE, indeksy BLAST (unite_blast_db) oraz Kraken2
├── dane_surowe/                # Wejściowe pliki .fastq.gz (np. barcode25.fastq.gz)
├── dane_posrednie/             # Folder tymczasowy na klastry i pliki pośrednie
├── skrypty/                    # Skrypty wykonawcze i analityczne
│   ├── blast_app.py            # Aplikacja Streamlit do wizualizacji wyników
│   ├── master_pipeline.sh      # Wersja testowa potoku
│   ├── process_all.sh          # Główny potok przetwarzania danych dla wszystkich próbek
│   ├── run_blast.sh            # Skrypt do lokalnej identyfikacji BLAST
│   ├── unite_to_kraken.py      # Konwerter bazy UNITE do formatu Kraken2
│   └── update_unite.sh         # Skrypt budowy/aktualizacji bazy Kraken2
├── wyniki_blast/               # Wyniki identyfikacji w formacie .tsv
├── wyniki_kraken/              # Raporty klasyfikacji Kraken2 (podzielone na próbki)
└── wyniki_spoa_konsensus/      # Ostateczne sekwencje FASTA dla każdego klastra
```

## Instrukcja Obsługi

### 1. Budowa środowiska
Zbuduj kontener Docker (wykonaj raz lub po jakiejkolwiek zmianie w pliku Dockerfile):
```bash
docker build -t grzyby_pro .
```

### 2. Przygotowanie bazy Kraken2 (Tylko za pierwszym razem)
Zanim potok użyje Krakena2 do wstępnej klasyfikacji, baza UNITE musi zostać przekonwertowana i zbudowana. Uruchom dedykowany skrypt przygotowujący bazę:
```bash
docker run --rm -v ${PWD}:/data grzyby_pro bash /data/skrypty/update_unite.sh
```

### 3. Przygotowanie lokalnej bazy BLAST (Tylko za pierwszym razem)
Zbuduj indeksy z pliku FASTA bazy UNITE dla precyzyjnego przyrównania końcowego:
```bash
docker run --rm -v ${PWD}:/data grzyby_pro makeblastdb \
  -in /data/baza_ITS/sh_general_release_dynamic_19.02.2025.fasta \
  -dbtype nucl -out /data/baza_ITS/unite_blast_db
```

### 4. Przetwarzanie surowych danych
Uruchom główny potok dla wszystkich próbek w folderze `dane_surowe/`. Skrypt wykona wycinanie adapterów (Porechop), filtrację długości (>300bp), klastrowanie (CD-HIT na poziomie 98%), złoży sekwencje konsensusowe (SPOA) i dokona wstępnej taksonomii (Kraken2):
```bash
docker run --rm -v ${PWD}:/data grzyby_pro bash /data/skrypty/process_all.sh
```

### 5. Identyfikacja BLAST
Uruchom dopasowanie wszystkich wygenerowanych konsensusów do bazy UNITE. Skrypt wygeneruje gotowe tabele `.tsv` w folderze `wyniki_blast/`:
```bash
docker run --rm -v ${PWD}:/data grzyby_pro bash /data/skrypty/run_blast.sh
```

### 6. Uruchomienie Dashboardu (Wizualizacja analityczna)
Aby odpalić aplikację webową, która automatycznie przetwarza tabele i rysuje statystyki, uruchom kontener z otwartym portem sieciowym (8501):
```bash
docker run -it --rm -v ${PWD}:/data -p 8501:8501 grzyby_pro streamlit run /data/skrypty/blast_app.py
```
Aplikacja będzie dostępna w przeglądarce internetowej pod adresem: **http://localhost:8501**

## Kluczowe Metryki w Dashboardzie
* **Percent Identity (pident):** Wskazuje na pewność taksonomiczną (przyjęto gatunek >97%, rodzaj >90%).
* **E-value:** Statystyczna istotność dopasowania (im niższa wartość, zbliżona do zera, tym pewniejszy wynik).
* **Top 1 Hit:** Wyczyszczona i podzielona systematyka z bazy UNITE (od Królestwa do Gatunku) dla najlepszego dopasowania każdego klastra.

---
