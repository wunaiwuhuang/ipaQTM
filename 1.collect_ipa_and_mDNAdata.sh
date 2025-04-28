#!/bin/bash

# 定义组织名称列表
cancers=("BLCA" "BRCA" "COAD" "DLBC" "ESCA" "GBM" "HNSC" "KICH" "KIRC" "KIRP" "LAML" "LGG" "LIHC" "LUAD" "LUSC" "MESO" "OV" "PAAD" "PCPG" "PRAD" "READ" "SARC" "SKCM" "STAD" "TGCT" "THCA" "THYM" "UCEC" "UCS" "UVM" "ACC" "CESC" "CHOL")


# 定义源目录和目标目录
src_base1="/data1/liuxiaochuan/project/TCGA_IPA2/01.data_prepare2/04.TCGA_usage"
src_base2="/data1/liuxiaochuan/project/TWAS_IPA/01.QTL/TCGA_hg38"
src_base3="/data1/wangwenhui/pu-DNAm/TCGA-GDC-DNAm"
dest_base="/data1/wuguojia/data/IPA_QTM_tcga/data"

# 确保目标目录存在
mkdir -p "$dest_base"

# ipa usage
    for cancer in "${cancers[@]}"; do
        src_dir="$src_base1/$cancer"
        # 确保源目录存在
        if [ -d "$src_dir" ]; then
            for file in "${cancer}.IPA_usage.txt"; do
                src_file="$src_dir/$file"
                dest_file="$dest_base/${cancer}_ipause.txt"
                
                # 检查文件是否存在，然后复制
                if [ -f "$src_file" ]; then
                    cp "$src_file" "$dest_file"
                    echo "Copied: $src_file -> $dest_file"
                else
                    echo "Warning: $src_file does not exist."
                fi
            done
        else
            echo "Warning: Directory $src_dir does not exist."
        fi
    done

# ipa location
    for cancer in "${cancers[@]}"; do
        src_dir="$src_base2/$cancer/Matrix_iQTL"
        # 确保源目录存在
        if [ -d "$src_dir" ]; then
            for file in "IPA_location.txt"; do
                src_file="$src_dir/$file"
                dest_file="$dest_base/${cancer}_ipaloc.txt"
                
                # 检查文件是否存在，然后复制
                if [ -f "$src_file" ]; then
                    cp "$src_file" "$dest_file"
                    echo "Copied: $src_file -> $dest_file"
                else
                    echo "Warning: $src_file does not exist."
                fi
            done
        else
            echo "Warning: Directory $src_dir does not exist."
        fi
    done

#mdna usage
    for cancer in "${cancers[@]}"; do
        src_dir="$src_base3/TCGA-${cancer}"
        # 确保源目录存在
        if [ -d "$src_dir" ]; then
            for file in "TCGA-${cancer}.DNAm.primary.uniq.matrix"; do
                src_file="$src_dir/$file"
                dest_file="$dest_base/${cancer}_mdnause.txt"
                
                # 检查文件是否存在，然后复制
                if [ -f "$src_file" ]; then
                    cp "$src_file" "$dest_file"
                    echo "Copied: $src_file -> $dest_file"
                else
                    echo "Warning: $src_file does not exist."
                fi
            done
        else
            echo "Warning: Directory $src_dir does not exist."
        fi
    done

echo "All files processed."

