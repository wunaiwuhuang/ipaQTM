# coding: utf-8
# i download protein annotation from uniprot https://www.uniprot.org/help/downloads, the uniprotKB-reviewed(swiss-prot) file. this script aims to extract useful infomation from uniprot_sprot.dat
input_file = "uniprot_sprot.dat"
id_file = "protein_coding_ensemblid.txt"
output_file = "uniprot_sprot_foruse.txt"
#load ensembl set
with open(id_file) as f:
    id_set = set(line.strip() for line in f if line.strip())
#deal with protein annotation file
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    entry_lines = []
    for line in infile:
        entry_lines.append(line)
        if line.startswith("//"):
            entry_text = ''.join(entry_lines)
            if any(ensembl_id in entry_text for ensembl_id in id_set):
                outfile.write(entry_text)
            entry_lines = []