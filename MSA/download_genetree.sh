#!/bin/bash
# MSA are downloaded in the folder of genetree, please replace the folder name as you wish

if [ "$#" -ne 1 ]; then
	input="ensembl_gene_id.txt"
else
	input=$1
fi

while IFS=$'\t' read -r gene name
do
	curl 'http://rest.ensembl.org/genetree/member/id/'"$gene"'?aligned=1' -H 'Content-type:application/json' >genetree/$gene.json

	result=`grep 'error\|problem' genetree/$gene.json`

	if [ -n "$result" ]
	then
		rm genetree/$gene.json
		echo $gene
	fi

done < "$input"

