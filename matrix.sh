#!/bin/bash

# 指定要合并的文件列表
files=()
for i in $(seq 58 81); do
    files+=("SRR21763${i}.count.txt")
done

# 创建一个新文件来存储合并后的结果
output_file="merged_matrix.txt"
csv_output_file="merged_matrix.csv"

# 处理第一个文件，提取第一列作为行表头，并将其内容复制到新文件中
tail -n +2 "${files[0]}" | cut -f 1 > gene_ids.txt

# 添加文件名作为标题到新文件
echo -e "geneid\t$(IFS=$'\t'; echo "${files[*]}")" > "$output_file"

# 创建临时文件存储列数据
tmp_output="tmp_output.txt"
paste gene_ids.txt > "$tmp_output"

# 处理每个文件，将各个样本的数据添加到新文件中
for file in "${files[@]}"; do
    paste "$tmp_output" <(tail -n +2 "$file" | cut -f 2) > temp.txt
    mv temp.txt "$tmp_output"
done

# 最终输出结果合并到output_file
cat "$tmp_output" >> "$output_file"
rm gene_ids.txt "$tmp_output"

# 移除第一行（文件名标题行）
tail -n +2 "$output_file" > temp.txt && mv temp.txt "$output_file"

# 将merged_matrix.txt转换为CSV文件
tr '\t' ',' < "$output_file" > "$csv_output_file"


