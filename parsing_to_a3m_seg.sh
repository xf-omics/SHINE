#!/bin/bash

if [ "$#" -ne 1 ]; then
	input="ensembl_gene_id.txt"
else
	input=$1
fi
			
while IFS= read -r gene
do
	file_input="../a3m_trim/"$gene.a3m
	if [ -f "$file_input" ];
	then
		./parsing_to_a3m_seg.o $gene
	fi
done < "$input"

