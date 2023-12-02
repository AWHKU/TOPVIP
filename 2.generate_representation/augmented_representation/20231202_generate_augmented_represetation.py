#conda activate py39
import numpy as np
import pandas as pd
import pickle
#### MLDE run on python 3.6.8 so make sure picke save at 4
from Bio import SeqIO


ALL_AAS = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
           "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-")
ALLOWED_AAS = set(ALL_AAS)
one_hot_dict = {aa: i for i, aa in enumerate(ALL_AAS)}
def convert_aa_to_one_hot(aa):
    onehot_array=np.zeros([1, 21])
    onehot_array[0, aa]=1
    return onehot_array

def load_fasta(data_filename):
    """ generate a list of one-hot array seq """
    """generate ids and labels list """
    ids = []
    labels = []
    seqs = []
    for record in SeqIO.parse(data_filename, "fasta"):
        ids.append(record.id)
        labels.append(record.description)
        a=np.array([convert_aa_to_one_hot(one_hot_dict[aa]) for aa in str(record.seq).upper().replace('.', '-')])
        seqs.append(a.reshape(-1))
    seqs=np.array(seqs)
    return ids, labels, seqs

def generate_augmented_representation_wrapper(fasta_filename, vae_filename, comb_ind_file, new_vae_encoding, new_comb_file, EV_features_files, new_EV_encoding, spec_column="AACombo"):
    ids, labels, seqs=load_fasta(fasta_filename)
    features=np.load(vae_filename)
    with open(comb_ind_file, "rb") as handle:
        comb_ind=pickle.load(handle)
    arranged_feat=[features[comb_ind[aa], ] for aa in ids]
    arranged_feat=np.array(arranged_feat)
    new_vae=np.column_stack((seqs, arranged_feat))
    new_comb_ind={k:v for v, k  in enumerate(ids)}
    np.save(new_vae_encoding, new_vae)
    with open(new_comb_file, 'wb') as handle:
        pickle.dump(new_comb_ind, handle, protocol=pickle.DEFAULT_PROTOCOL)
    ### rearrange vae feautures based on ids 
    EV_features=pd.read_csv(EV_features_files)
    EV_features=EV_features.set_index(spec_column, drop=False)
    new_plmc=EV_features.reindex(ids)
    new_plmc=new_plmc["plmc"].values
    new_EV=np.column_stack((seqs, new_plmc))
    np.save(new_EV_encoding, new_EV)
    return seqs, new_EV, new_vae, new_comb_ind


### SpCas9 
data_filename="data/generate_representations/SpCas9_augmented/SpCas9_AACombo_only.fasta"
vae_filename="data/generate_representations/SpCas9_VAE/SpCas9variants_AACombo_VAE_embeddedZ30.npy"
comb_ind_file="data/generate_representations/SpCas9_VAE/SpCas9variants_AACombo_VAEz30_ComboToIndex.pkl"
new_vae_encoding="SpCas9_augmentedVAE_embedded.npy"
new_comb_file="SpCas9_augmentedVAE_ComboToIndex.pkl"
EV_features_files="data/generate_representations/SpCas9_augmented/data_EV.csv"
new_EV_encoding="SpCas9_augmentedEVmutation_embedded.npy"
seqs, new_EV_encoding, new_vae_encoding, new_comb_ind=generate_augmented_representation_wrapper(data_filename, vae_filename, comb_ind_file, new_vae_encoding, new_comb_file, EV_features_files, new_EV_encoding, spec_column="AACombo")

### SaCas9 1296 
data_filename="data/generate_representations/SaCas9_augmented/887888887985986988989991_AACombo_only.fasta"
vae_filename="data/generate_representations/SaCas9_VAE/SaCas9_variants_AACombo_VAEz30_embedded.npy"
comb_ind_file="data/generate_representations/SaCas9_VAE/SaCas9_variants_AACombo_VAEz30_ComboToIndex.pkl"
new_vae_encoding="SaCas9_1296_augmentedVAE_embedded.npy"
new_comb_file="SaCas9_1296_augmentedVAE_ComboToIndex.pkl"
EV_features_files="data/generate_representations/887888887985986988989991_data_EV.csv"
new_EV_encoding="SaCas9_1296_augmentedEVmutation_embedded.npy"
seqs, new_EV_encoding, new_vae_encoding, new_comb_ind=generate_augmented_representation_wrapper(data_filename, vae_filename, comb_ind_file, new_vae_encoding, new_comb_file, EV_features_files, new_EV_encoding, spec_column="AACombo")

### SaCas9 8000
data_filename="data/generate_representations/SaCas9_augmented/888889909_AACombo_only.fasta"
vae_filename="data/generate_representations/SaCas9_VAE/SaCas9_888889909_AACombo_VAEz30_embedded.npy"
comb_ind_file="data/generate_representations/SaCas9_VAE/SaCas9_888889909_AACombo_VAEz30_ComboToIndex.pkl"
new_vae_encoding="SaCas9_888889909_augmentedVAE_embedded.npy"
new_comb_file="SaCas9_888889909_augmentedVAE_ComboToIndex.pkl"
EV_features_files="data/generate_representations/SaCas9_augmented/888889909_data_EV.csv"
new_EV_encoding="SaCas9_888889909_augmentedEVmutation_embedded.npy"
seqs, new_EV_encoding, new_vae_encoding, new_comb_ind=generate_augmented_representation_wrapper(data_filename, vae_filename, comb_ind_file, new_vae_encoding, new_comb_file, EV_features_files, new_EV_encoding, spec_column="AACombo")


### PhoQ -- 140517 entries from PhoQ.xlsx used in Clade 2.0 
fasta_filename="data/generate_representations/PhoQ_augmented/PhoQ_aacombo_only.fasta"
vae_filename="data/generate_representations/PhoQ_VAE/PhoQ_variants_AACombo_VAEz30_embedded.npy"
comb_ind_file="data/generate_representations/PhoQ_VAE/PhoQ_variants_AACombo_VAEz30_ComboToIndex.pkl"
new_vae_encoding=/PhoQ_augmentedVAE_embedded.npy"
new_comb_file="PhoQ_augmentedVAE_ComboToIndex.pkl"
EV_features_files="data/generate_representations/PhoQ_augmented/PhoQ_fitness_EV.csv"
new_EV_encoding="PhoQ_augmentedEVmutation_embedded.npy"
seqs, new_EV_encoding, new_vae_encoding, new_comb_ind=generate_augmented_representation_wrapper(fasta_filename, vae_filename, comb_ind_file, new_vae_encoding, new_comb_file, EV_features_files, new_EV_encoding, spec_column="Variants")

### PABP
fasta_filename="data/generate_representations/PABP_augmented/PABP_11-87-mutations.fasta"
vae_filename="data/generate_representations/PABP_VAE/PABP_variants_AACombo_VAEz30_embedded.npy"
comb_ind_file="data/generate_representations/PABP_VAE/PABP_variants_AACombo_VAEz30_ComboToIndex.pkl"
new_vae_encoding="PABP_augmentedVAE_embedded.npy"
new_comb_file="PABP_augmentedVAE_ComboToIndex.pkl"
EV_features_files="data/generate_representations/PABP_fitness_EV_rename.csv"
new_EV_encoding="PABP_augmentedEVmutation_embedded.npy"
seqs, new_EV_encoding, new_vae_encoding, new_comb_ind=generate_augmented_representation_wrapper(fasta_filename, vae_filename, comb_ind_file, new_vae_encoding, new_comb_file, EV_features_files, new_EV_encoding, spec_column="mutations")

### UBE 
fasta_filename="data/generate_representations/UBE_augmented/UBE4B/UBE_23-100-mutations.fasta"
vae_filename="data/generate_representations/UBE_VAE/UBE_variants_AACombo_VAEz30_embedded.npy"
comb_ind_file="data/generate_representations/UBE_VAE/UBE_variants_AACombo_VAEz30_ComboToIndex.pkl"
new_vae_encoding="UBE_augmentedVAE_embedded.npy"
new_comb_file="UBE_augmentedVAE_ComboToIndex.pkl"
EV_features_files="data/generate_representations/UBE_augmented/UBE_fitness_EV_rename.csv"
new_EV_encoding="UBE_augmentedEVmutation_embedded.npy"
seqs, new_EV_encoding, new_vae_encoding, new_comb_ind=generate_augmented_representation_wrapper(fasta_filename, vae_filename, comb_ind_file, new_vae_encoding, new_comb_file, EV_features_files, new_EV_encoding, spec_column="mutations")

### GFP 
fasta_filename="data/generate_representations/GFP_augmented/GFP_14-151-mutations.fasta"
vae_filename="data/generate_representations/GFP_VAE/GFP_variants_AACombo_VAEz30_embedded.npy"
comb_ind_file="data/generate_representations/GFP_VAE/GFP_variants_AACombo_VAEz30_ComboToIndex.pkl"
new_vae_encoding="GFP_augmentedVAE_embedded.npy"
new_comb_file="GFP_augmentedVAE_ComboToIndex.pkl"
EV_features_files="data/generate_representations/GFP_augmented/GFP_fitness_EV_rename.csv"
new_EV_encoding="GFP_augmentedEVmutation_embedded.npy"
seqs, new_EV_encoding, new_vae_encoding, new_comb_ind=generate_augmented_representation_wrapper(fasta_filename, vae_filename, comb_ind_file, new_vae_encoding, new_comb_file, EV_features_files, new_EV_encoding, spec_column="mutations")
