#!/bin/bash

if [ "$#" -ne 1 ]; then
	input="../ensembl_gene_id.txt"
else
	input=$1
fi
			
while IFS= read -r gene
do
	file_input="../msa/"$gene.txt
	if [ -f "$file_input" ];
	then
		./parsing_to_a3m.o $gene
		cat temp >> ../a3m/$gene.a3m
		rm temp
	fi
done < "$input"

