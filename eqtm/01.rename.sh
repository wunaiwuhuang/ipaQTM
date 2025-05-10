#!/bin/bash

# cp -r /data1/wangwenhui/pu-DNAm/eQTM/01.input/DNAm_location  ./
# cp -r /data1/wangwenhui/pu-DNAm/eQTM/01.input/DNAm_matrix ./
# cp -r /data1/wangwenhui/pu-DNAm/eQTM/01.input/TPM_location    ./
# cp -r /data1/wangwenhui/pu-DNAm/eQTM/01.input/TPM_matrix  ./
# cp -r /data1/wangwenhui/pu-DNAm/eQTM/01.input/Covariate_matrix    ./

cd /data1/wuguojia/data/IPA_QTM_tcga/eqtm/data
mv ./Covariate_matrix/* ./
mv ./DNAm_location/* ./
mv ./DNAm_matrix/* ./
mv ./TPM_location/* ./
mv ./TPM_matrix/* ./

rm -rf ./Covariate_matrix ./DNAm_location ./DNAm_matrix ./TPM_location ./TPM_matrix

for file in *.cov.matrix; do
  newname="${file/.cov.matrix/_covariates.txt}"
  mv "$file" "$newname"
done

for file in *.DNAm.location.txt; do
  newname="${file/.DNAm.location.txt/_mdnaloc.txt}"
  mv "$file" "$newname"
done

for file in *.DNAm.matrix; do
  newname="${file/.DNAm.matrix/_mdnause.txt}"
  mv "$file" "$newname"
done

for file in *.tpm.location.txt; do
  newname="${file/.tpm.location.txt/_tpmloc.txt}"
  mv "$file" "$newname"
done

for file in *.tpm.matrix; do
  newname="${file/.tpm.matrix/_tpmuse.txt}"
  mv "$file" "$newname"
done