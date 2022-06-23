#!/bin/bash

if [ "$#" -ne 1 ]; then
	input="ensembl_genetree_gene_protein_id.txt"
else
	input=$1
fi
			
th=300

while IFS=$'\t' read -r -a array
do
	gene=${array[0]}
	file_input="../a3m/"$gene.a3m
	if [ "${array[2]}" -gt ${th} ];
	then
		./parsing_to_short_a3m.o $gene $th
		./remove_gaps.o $gene
		cp temp ../a3m_trim/$gene.a3m
	else
		cp $file_input ../a3m_trim/
	fi
done < "$input"

