#!/bin/bash
cd ./data/
cancers=()
for file in *_ipause.txt; do
    [[ -e "$file" ]] || continue
    cancer_name="${file%%_ipause.txt}" #remove suffix
    cancers+=("$cancer_name")
done

if [[ ${#cancers[@]} -eq 0 ]]; then
    echo "未找到匹配的文件，退出。"
    exit 1
fi

read -p "请输入 threshold 值: " threshold #input threshold

filtered_cancers=()
for cancer in "${cancers[@]}"; do
    ipa_file="${cancer}_ipause.txt"
    if [[ ! -f "$ipa_file" ]]; then
        echo "警告: $ipa_file 不存在，跳过 $cancer"
        continue
    fi
    num_columns=$(head -n 1 "$ipa_file" | awk -F'\t' '{print NF}') #calculate the number of columns

    if (( num_columns >= threshold + 1 )); then
        filtered_cancers+=("$cancer")
    else
        echo "移除 $cancer (列数: $num_columns，小于 ${threshold}+1)"
    fi
done

removed_count=$(( ${#cancers[@]} - ${#filtered_cancers[@]} )) #calculate remove cancer count
remaining_count=${#filtered_cancers[@]} #calculate remaining cancer count

echo "最终筛选后的癌症列表（剩余数量: $remaining_count，抹去数量: $removed_count）:"
for cancer in "${filtered_cancers[@]}"; do
    ipa_file="${cancer}_ipause.txt"
    num_columns=$(head -n 1 "$ipa_file" | awk -F'\t' '{print NF}') #calculate the number of columns
    echo "$cancer (列数: $num_columns)"
done

#delect the cancer file
# ask for confirmation before deleting
read -p "是否删除不符合条件的癌症文件？(y/n): " confirm
if [[ "$confirm" != "y" ]]; then
    echo "未删除文件。"
    exit 0
fi
# delete files for cancers not in filtered_cancers
for cancer in "${cancers[@]}"; do
    if [[ ! " ${filtered_cancers[@]} " =~ " ${cancer} " ]]; then
        rm -rf "${cancer}_ipause.txt"
        rm -rf "${cancer}_ipaloc.txt"
        rm -rf "${cancer}_mdnause.txt"
        echo "删除 $cancer 的文件"
    fi
done