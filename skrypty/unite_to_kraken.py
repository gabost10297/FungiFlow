import sys
import os

if len(sys.argv) != 2:
    print("BŁAD: Podaj plik wejściowy.")
    sys.exit(1)

input_fasta = sys.argv[1]
out_fasta = "/data/baza_ITS/unite_kraken.fasta"
out_dir = "/data/baza_ITS/taxonomy"
os.makedirs(out_dir, exist_ok=True)

print("Rozpoczynam tłumaczenie bazy UNITE na język Kraken2...")

names_file = open(os.path.join(out_dir, "names.dmp"), "w")
nodes_file = open(os.path.join(out_dir, "nodes.dmp"), "w")
fasta_file = open(out_fasta, "w")

# Root taksonomii
nodes_file.write("1\t|\t1\t|\tno rank\t|\t\n")
names_file.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")

tax_dict = {"root": 1}
current_id = 2
counter = 1 # Licznik dla unikalności ID

rank_map = {'k': 'kingdom', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}

with open(input_fasta, "r") as f:
    for line in f:
        if line.startswith(">"):
            header = line.strip()[1:]
            # Wyciągamy pierwszy człon nazwy grzyba
            clean_name = header.replace('|', ' ').split()[0]
            # Tworzymy unikalne ID dla Krakena
            seq_id = f"{clean_name}_{counter}"
            counter += 1
            
            parts = header.split('|')
            tax_str = next((p for p in parts if "k__" in p), "")
            parent_id = 1 
            
            if tax_str:
                lineages = tax_str.split(';')
                for tax in lineages:
                    if not tax or "__" not in tax: continue
                    res = tax.split('__', 1)
                    if len(res) < 2: continue
                    rank_code, name = res[0], res[1]
                    if name.lower() == "unidentified" or not name: continue
                    
                    node_path = f"{parent_id}_{name}"
                    if node_path not in tax_dict:
                        tax_dict[node_path] = current_id
                        rank_name = rank_map.get(rank_code, "no rank")
                        nodes_file.write(f"{current_id}\t|\t{parent_id}\t|\t{rank_name}\t|\t\n")
                        names_file.write(f"{current_id}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
                        parent_id = current_id
                        current_id += 1
                    else:
                        parent_id = tax_dict[node_path]
            
            fasta_file.write(f">{seq_id}|kraken:taxid|{parent_id}\n")
        else:
            fasta_file.write(line)

names_file.close()
nodes_file.close()
fasta_file.close()
print(f"SUKCES! Przygotowano {counter-1} sekwencji.")