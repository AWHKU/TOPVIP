# TopVIP
Scripts for running the analyses in the paper **"Accurate top protein variant discovery via low-N pick-and-validate machine learning"**.
There are in total 7 parts of the analyses, from generating zero-shot predictions to final plotting. 
The input, intermediate, and output data could be found in Figshare , DOI: https://doi.org/10.6084/m9.figshare.24715761

## 1. Generating Zero-shot prediction 
The first step of the analysis was to generate zero-shot scores of the protein variants. We generated zero-shot prediction scores from Prot-Bert, MSA "esm_msa1_t12_100M_UR50S", EVmutation, Seqdesign, arDCA, EVoEF DDG, and ESM2 scores and voting from efficient evolution method. Go into the subfolders of *1.generate_zeroshot_predictors* for details to generate each type of prediction scores. 

### Prot-bert and MSA transformer models
We use predict_zero_shot.py from [MLDE](https://github.com/fhalab/MLDE.git) to generate zero-shot prediction scores from pre-built model "esm_msa1_t12_100M_UR50S" and "prot_bert". 

Follow the instruction in MLDE github section **"Zero Shot Prediction with predict_zero_shot.py"** to generate the scores.
#### Reference 
Wittmann et al. Informed training set design enables efficient machine learning-assisted directed protein evolution. DOI:https://doi.org/10.1016/j.cels.2021.07.008
### EVmutation

We used the [EVcoupling webserver](https://v2.evcouplings.org/) to generate the MSA and coupling models for SpCas9 and SaCas9. We used the environment and supporting src in git hub repository [combining-evolutionary-and-assay-labelled-data](https://github.com/chloechsu/combining-evolutionary-and-assay-labelled-data.git), follow the instruction there for setting up the environment to run the scripts. 
The benchmarking datasets of GFP, PABP, and UBE are derived from Hsu et al [https://doi.org/10.6078/D1K71B]( https://doi.org/10.6078/D1K71B). 
The benchmarking dataset PhoQ are derived from [CLADE2.0](https://github.com/YuchiQiu/CLADE-2.0.git); use the download_data.sh script to download the dataset. 
#### Reference
Hopf et al. Mutation effects predicted from sequence co-variation. DOI: https://doi.org/10.1038/nbt.3769

Hopf et al. The EVcouplings Python framework for coevolutionary sequence analysis. DOI: https://doi.org/10.1093/bioinformatics/bty862

Hsu et al. Learning protein fitness models from evolutionary and assay-labeled data. DOI: https://doi.org/10.1038/s41587-021-01146-5

Qiu and Wei. CLADE 2.0: Evolution-Driven Cluster Learning-Assisted Directed Evolution. DOI: https://doi.org/10.1021/acs.jcim.2c01046

### Seqdesign

We run Seqdesign based on the default setting provided in [seqdesign-pytorch github](https://github.com/aaronkollasch/seqdesign-pytorch.git). 
#### Reference
Shin et al. Protein design and variant prediction using autoregressive generative models, DOI: https://doi.org/10.1038/s41467-021-22732-w
### arDCA

We followed the instruction provided in the [arDCA github](https://github.com/pagnani/ArDCA.jl.git) to run arDCA. The MSA alignments are generated or obtained from running EVcoupling described in previous section EVmutation. 
#### Reference
Trinquier et al. Efficient generative modeling of protein sequences using simple autoregressive models. DOI: https://doi.org/10.1038/s41467-021-25756-4
### EvoEF

We follow the instruction provided in the [EvoEF github](https://github.com/tommyhuangthu/EvoEF2.git)
 to run EvoEF. 
#### Reference
Huang et al. EvoEF2: accurate and fast energy function for computational protein design. DOI:https://doi.org/10.1093/bioinformatics/btz740 

### efficient evolution - ESM2 

We used the enviroment and scripts provided in the [efficient-evolution github](https://github.com/brianhie/efficient-evolution.git) for generating the ESM2 scores. 
#### Reference
Hie et al. Efficient evolution of human antibodies from general protein language models. DOI: https://doi.org/10.1038/s41587-023-01763-2
### Dataset Reference
Podgornaia and Laub. Pervasive degeneracy and epistasis in a protein-protein interface. DOI: https://doi.org/10.1126/science.1257360

Sarkisyan et al. Local fitness landscape of the green fluorescent protein. DOI: https://doi.org/10.1038/nature17995

Starita et al. Activity-enhancing mutations in an E3 ubiquitin ligase identified by high-throughput mutagenesis. DOI:https://doi.org/10.1073/pnas.1303309110

Melamed et al. Deep mutational scanning of an RRM domain of the Saccharomyces cerevisiae poly (A)-binding protein. DOI: https://doi.org/10.1261/rna.040709.113 

Thean et al. Machine learning-coupled combinatorial mutagenesis enables resource-efficient engineering of CRISPR-Cas9 genome editor activities. DOI: https://doi.org/10.1038/s41467-022-29874-5

Choi et al. Combinatorial mutagenesis en masse optimizes the genome editing activities of SpCas9 DOI: https://doi.org/10.1038/s41592-019-0473-0. 
 
## 2. Generate protein representation
We generated Georgiev representation based on the instruction of [MLDE](https://github.com/fhalab/MLDE.git).


We followed the instruction in 
distances.ipynb in the github repository [meaningful-protein-representations](https://github.com/MachineLearningLifeScience/meaningful-protein-representations.git) to generate the VAE latent representation. 


We generated the augmented representation from the Georgiev representation and concatenate the EVmutation score and the VAE representation. 
#### Reference
Detlefsen, Hauberg and Boomsma. Learning meaningful representations of protein sequences. DOI: https://doi.org/10.1038/s41467-022-29443-w

Hsu et al. Learning protein fitness models from evolutionary and assay-labeled data. DOI: https://doi.org/10.1038/s41587-021-01146-5

## 3. Clustered sampling 
We gnerated K-means clustering based on the Georgiev representation of the protein. We adopt the clustering_sampling.py from [CLADE github](https://github.com/YuchiQiu/CLADE.git) to perform clustering, ranking and adaptive sampling for the preparation of MLDE input Fitness files. 
#### Reference
Qiu, Hu and Wei. Cluster learning-assisted directed evolution. DOI: https://doi.org/10.1038/s43588-021-00168-y

## 4. Multi-round sampling (TopVIP)
The R scripts use the zero-shot predicted scores or MLDE results of previous rounds to prepare the input MLDE Fitness files for the next round of MLDE. 

We run TopVIP on both in-house and benchmarking datasets.

## 5. Running CLADE 2.0
We run CLADE 2.0 using the zero-shot predictor scores only in the benchmarking section. We modified the clustering_sampling.py in [CLADE 2.0](https://github.com/YuchiQiu/CLADE-2.0.git) to remove the codes that were required in the empirical data of the full library as the input. 

#### Reference
Qiu and Wei. CLADE 2.0: Evolution-Driven Cluster Learning-Assisted Directed Evolution. DOI: https://doi.org/10.1021/acs.jcim.2c01046

## 6. Summarize MLDE performance
The R scripts were used to evaluate the performance of the MLDE runs. 

## 7. Plotting
The R scripts were used to plot the figures and supplementary figures in the paper.
#### Reference
R Core Team (2021). R: A language and environment for statistical computing. https://www.R-project.org/

Wickham et al. Welcome to the tidyverse. DOI: https://doi.org/10.21105/joss.01686

Wickham et al. ggplot2: Elegant Graphics for Data Analysis. https://ggplot2.tidyverse.org

Wickham and Byran. readxl: Read Excel Files. https://CRAN.R-project.org/package=readxl

Slowikowski et al. ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. https://CRAN.R-project.org/package=ggrepel

