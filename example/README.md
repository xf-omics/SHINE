# SHINE examples
A computational method to predict pathogenicity for inframe indels.

SHINE leverages pre-trained protein language models and limited available pathogenicity labels for infram indels to improve variant interpretation using transfer learning.

## 1. Variants
run SHINE on variants with known pathogenicity labels
bash run_variants.sh


## 2. Gene
run SHINE on an entire gene without known pathogenicity labels
bash run_gene.sh

### trouble shooting
If the PCA eigenvectors cannot be loaded, please use scilit-learn version 1.1.1.
