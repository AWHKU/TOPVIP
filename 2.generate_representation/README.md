# Generate VAE representation from Proteins 

## VAE
We used the multiple sequence alignments (MSA) of the benchmarking proteins:PhoQ, PABP, GFP, and UBE from [CLADE2.0](https://github.com/YuchiQiu/CLADE-2.0.git) and 
Hsu et al [chloechsu/combining-evolutionary-and-assay-labelled-data](https://github.com/chloechsu/combining-evolutionary-and-assay-labelled-data.git).

For the in-house libraries of SpCas9 and SaCAs9, the MSA input could be found in data/zeroshot_prediction/EVcoupling. 

To generate the VAE representation of the protein, we adopt the distances.ipynb from [https://github.com/MachineLearningLifeScience/meaningful-protein-representations.git](https://github.com/MachineLearningLifeScience/meaningful-protein-representations.git).


Run each of the python script in ./VAE to generate the relevant VAE representation. The output are stored at data/generate_representations. 

## Augmented representation 
We built up upon the Georgiev representation generated from MLDE [generate_encoding.py](https://github.com/fhalab/MLDE/blob/main/generate_encoding.py) and combine the EVmutation score and the VAE representation with the Georgiev representation. 


The input of the python script are stored at data/generate_representations.
Run the python script 20231202_generate_augmented_represetation.py to generate the augmented representations. 
