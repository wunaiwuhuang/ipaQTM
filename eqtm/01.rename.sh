#!/bin/bash

cd /data1/wuguojia/data/IPA_QTM_tcga/eqtm/data
cp -r /data1/wangwenhui/pu-DNAm/eQTM/TCGA-Xena-expression/*.tsv    ./

for file in *.tpm.tsv; do
  prefix=$(echo "$file" | sed -r 's/^TCGA-([^.]+)\.tpm\.tsv$/\1/')
  newname="${prefix}_tpmuse.txt"
  mv "$file" "$newname"
  awk 'NR==1 {$1="id"} 1' OFS="\t" "$newname" > tmp && mv tmp "$newname" # send first element to id
done