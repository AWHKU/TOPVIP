#Clustered sampling#
The python scripts were used for clustering the Georgiev representation of the protein tested in the paper. 

In brief, it used the representation for K-means clustering. The number of clusters were determined manually based on the elbow-plots.

The cluster allocation are used in subsequent sampling steps in Multi-round sampling and TOPVIP. 

Adaptive sampling were performed for SpCas9, SaCas9 WED-PI and SaCas9 WED libraries to evaluate the benefits of adaptive sampling with clustered-sampling.


The Georgiev representation were generated using MLDE generate_embedding.py (https://github.com/fhalab/MLDE.git). The output files *_georgiev_embedding.npy* and *_georgiev_ComboToIndex.pkl* are used as the input for the python scripts.


The input, intermeditate, and output files are present in data/clustered_sampling
