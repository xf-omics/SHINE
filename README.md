# SHINE
A computational method to predict pathogenicity for inframe indels.

SHINE leverages pre-trained protein language models and limited available pathogenicity labels for infram indels to improve variant interpretation using transfer learning.

# 1. Datasets
*All datasets are provided in the dataset folder

Training dataset includes one-amino-acid inframe indels from ClinVar and one-amino-acid inframe indels with allele frequency > 0.1% in gnomAD.

Validation dataset includes multiple-amino-acid inframe indels from ClinVar and multiple-amino-acid inframe indels with allele frequency > 0.1% in gnomAD.

NDD test dataset includes _de novo_ inframe indels in NDD genes from NDD cases and germline one-to-three-amino-acid inframe indels in the same NDD genes from UK biobank controls.

Cancer mutational hotspot test dataset includes somatic inframe indels in the hotspots from cancer cases and germline one-to-three-amino-acid inframe indels in the same cancer genes from UK biobank controls.

# 2. Pretraining
We used ESM-1b and MSA transformers to generate latent representations of amino acids for each protein

## 2.0 preprocessing for protein sequences
*Please note the maximum length of input proteins is 1024. I provide the script to divide protein sequences or MSA into pieces if length > 1024. Feel free to divide in the way you wish. Protein sequences are downloaded from ensembl and are not provided.

### ESM-1b deletion sequences
./esm1b_input_seq.o protein_transcript_list

### ESM-1b insertion sequences
This script generates protein sequences with inserted amino acids and change amino acid positions based on mutated protein sequences for the given variants at the same time.
./parsing_esm1b_input.o insertion_variant_list

### ESM-MSA sequences
This script is in the MSA directory.
./parsing_to_a3m_seg.o gene_accession_list

## 2.1 run ESM-1b model
python repr_esm1b.py esm1b_t33_650M_UR50S input_fasta output_dir --repr_layers 33 --include per_tok

## 2.2 run MSA transformer
python repr_msa_indel.py esm_msa1b_t12_100M_UR50S input_gene_accession_ids input_msa_dir output_dir --repr_layers 12

## 2.3 preprocessing for variants
### ESM-1b deletion variants
This script works for both deletions and insertions. 
./parsing_esm1b_input.o deletion_variant_list

### ESM-MSA variants
The input deletion list and insertion list are the same for ESM-1b and ESM-MSA.
./parsing_esm_input_a3m_del.o deletion_variant_list ensembl_genetree_trim_gene_protein_index.txt a3m_trim_seg/
./parsing_esm_input_a3m_ins.o insertion_variant_list ensembl_genetree_trim_gene_protein_index.txt a3m_trim_seg/

# 3. Running prediction model
python variant_prediction_esm1b_msa_del.py training_variants testing_variants prediction_result msa_representation_dir esm1b_representation_dir

python variant_prediction_esm1b_msa_ins.py training_variants testing_variants prediction_result msa_representation_dir esm1b_representation_dir

# 4. Predictions of deletions
We provide SHINE prediction scores for all possible single amino acid deletions in human genes:
https://www.dropbox.com/s/tc12faehbl6rv7g/SHINE_deletion_scores.tar.gz?dl=0

