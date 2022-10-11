# Description about how we processed the MSA data
The folder names are hardcoded, feel free to replace them as you wish

## 0. ensembl gene and protein accession ids are provided
ensembl_gene_id.txt
ensembl_genetree_gene_protein_id.txt

## 1. download gene tree from ensembl as \*.json in genetree/
download_genetree.sh

## 2. parsing \*.json to generate \*.txt in msa/ including species name, gene accession and protein sequence
python parsing_json.py ensembl_gene_id.txt genetree/ msa/

## 3. parsing \*.txt to \*.a3m in a3m/ as fasta format
bash parsing_to_a3m.sh ensembl_gene_id.txt

## 4. trim the tree (maximum 300 species) from a3m/\*.a3m to a3m_trim/\*.a3m
bash parsing_to_short_a3m.sh

## 5. remove gaps at the beginning and end of MSA. copy back to a3m_trim/\*.a3m
bash parsing_to_short_a3m.sh

## 6. divide MSA into segments if size>1023 (limited by ESM-MSA) from a3m_trim/\*.a3m to a3m_trim_seg/\*.a3m
bash parsing_to_a3m_seg.sh

## 7. evaluate the MSA quality using a3m/\*.a3m
MSA depth before and after the trim is provided as ensembl_genetree_gene_protein_depth_updated.txt
