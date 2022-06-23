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

1.1 run ESM-1b model
python repr_esm1b.py esm1b_t33_650M_UR50S input_fasta output_dir --repr_layers 33 --include per_tok

*Please note the maximum length of input proteins is 1024. I provide the script to divide protein sequences into pieces if length > 1024. Feel free to divide in the way you wish. Protein sequences are downloaded from ensembl and are not provided.

./esm1b_input_seq.cpp protein_transcript_list

1.2 run MSA transformer
python repr_msa_indel.py esm_msa1b_t12_100M_UR50S gene_accession_id input_msa_dir output_dir --repr_layers 12
