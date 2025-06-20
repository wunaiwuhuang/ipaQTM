# coding: utf-8
# extract keywords from all protein genes
input_file = "uniprot_sprot_foruse.txt"
id_file = "protein_coding_ensemblid.txt"
output_file = "uniprot_sprot_final.txt"
# read ensembl id
with open(id_file) as f:
    id_set = set(line.strip() for line in f if line.strip())
# deal with kw
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    entry_lines = []
    for line in infile:
        entry_lines.append(line)
        if line.startswith("//"):
            entry_text = ''.join(entry_lines)
            matched_ids = [ensembl_id for ensembl_id in id_set if ensembl_id in entry_text]
            if matched_ids:
                # extract kw lines
                kw_lines = [l.strip() for l in entry_lines if l.startswith("KW")]
                # remove KW prefix and end dot
                kw_clean = ' '.join(l[5:].rstrip('.').strip() for l in kw_lines)
                outfile.write(">" + ";".join(matched_ids) + "\n")
                outfile.write(kw_clean + "\n")
            entry_lines = []