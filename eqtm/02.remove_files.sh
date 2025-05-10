#!/bin/bash
path_ipa="/data1/wuguojia/data/IPA_QTM_tcga/data"

cancer_list=$(ls "$path_ipa"/*_covariates_peer.txt 2>/dev/null | sed 's#.*/##' | sed 's/_covariates_peer.txt//' | sort -u)

cd /data1/wuguojia/data/IPA_QTM_tcga/eqtm/data || exit 1
# delete files that do not start with any of the cancer types
for file in *; do
    keep=false
    for cancer in $cancer_list; do
        if [[ "$file" == ${cancer}_* ]]; then
            keep=true
            break
        fi
    done
    if [ "$keep" = false ]; then
        echo "Deleting $file"
        rm -f "$file"
    fi
done
