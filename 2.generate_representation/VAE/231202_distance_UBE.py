# cd /data/users/athchu/meaning_rep
# conda activate /data/users/athchu/.conda/envs/protein_rep
# export 'PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:1000'
# CUDA_LAUNCH_BLOCKING=1 python
import numpy as np
import pytorch_lightning as pl
import torch
import torch.nn.functional as F
import os
import re
from Bio import SeqIO
import torchplot
import sys
#sys.path.append("meaning_rep/meaningful-protein-representations/models")
from vae_geometric import VAE, get_hparams, aa1_to_index, aa1, numeric_curve_optimizer


### def get_data
data_filename="data/zeroshot_prediction/EVcoupling/UBE_homologues_w1.fasta"
torch.cuda.empty_cache()
ids = []
labels = []
seqs = []
for record in SeqIO.parse(data_filename, "fasta"):
        ids.append(record.id)
        labels.append(record.description)
        seqs.append(np.array([aa1_to_index[aa] for aa in str(record.seq).upper().replace('.', '-')]))

seqs = torch.from_numpy(np.vstack(seqs))
data=seqs
labels = np.array(labels)

weights = None
one_hot1 = F.one_hot(seqs[:len(seqs)//2].long()).bool()
one_hot2 = F.one_hot(seqs[len(seqs)//2:].long()).bool()
one_hot = torch.cat([one_hot1, one_hot2]).to('cuda' if torch.cuda.is_available() else 'cpu')
assert(len(seqs) == len(one_hot))
del one_hot1
del one_hot2
### any symbols not aa ==0 in one_hot , ="A"
one_hot[seqs>19] = 0
### original size: [7844, 263 , 23]
flat_one_hot = one_hot.flatten(1)
weights = {}
similarity_threshold=[0.2]
for st in similarity_threshold:
    weights[str(st)] = []


weight_batch_size = 10
flat_one_hot = flat_one_hot.float()
for i in range(seqs.size(0) // weight_batch_size + 1):
    x = flat_one_hot[i * weight_batch_size : (i + 1) * weight_batch_size]
    similarities = torch.mm(x, flat_one_hot.T)
    lengths = (seqs[i * weight_batch_size : (i + 1) * weight_batch_size] <=19).sum(1).unsqueeze(-1).to('cuda' if torch.cuda.is_available() else 'cpu')
    for st in similarity_threshold:
        w = 1.0 / (similarities / lengths).gt(st).sum(1).float()
        weights[str(st)].append(w)
        # plt.hist(similarities/lengths)


for st in similarity_threshold:
    weights[str(st)] = torch.cat(weights[str(st)])
    #neff = weights.sum()

### initiate weight, one sequence, one weight in each st
#weights["0.1"].size()
#torch.Size([999])
torch.save(weights, "UBE_VAE/UBE_homologues_AACombo_weightsz30.pkl")

# Train/retrieve VAE models
import time
### seems didn't really given me the model back to do encoding, what to do?
cmd_args = ['-gpu', '1' if torch.cuda.is_available() else '0', 
                '-kl_warmup_steps', str(1), 
                "-epochs", str(100), 
                "-zdim", str(30),
                "-train_fraction", str(0.8),
                "-val_fraction", str(0.2), 
                "-seed", str(123), 
                "-bs", str(16),
                "-lr", str(0.0001),
                "-iwae_bound", str(False),
                "-mask_out_gaps", str(True),                
                "-sparsity_prior", str(False),
                "-sparsity_prior_lambda", str(0.0001),
                "-simplify_to_ae", str(False)]

hparams = get_hparams(cmd_args)
pl.seed_everything(hparams.seed)
perm = np.random.permutation(data.shape[0])
print('Training model!')
model = VAE(data=data, weights=weights['0.2'], perm=perm, hparams=hparams)    


from pytorch_lightning.callbacks.early_stopping import EarlyStopping
trainer = pl.Trainer(accelerator='gpu', devices=1,
                             max_epochs=hparams.epochs,
                             callbacks=[EarlyStopping(monitor="val_loss", mode="min")])
### alternative early stopping
#early_stop_callback = EarlyStopping(monitor="val_acc", min_delta=0.00, patience=3, verbose=False, mode="max")
#trainer = Trainer(accelerator='gpu', devices=1,
#                             max_epochs=hparams.epochs,
#                             callbacks=[early_stop_callback])
#### training step ?
trainer.fit(model)

# test the model
#trainer.test(model, data=test_data?)
trainer.save_checkpoint("UBE_VAE/UBE_homologues_modelz30.ckpt")
encoded_train=model.embedding(data)
encoded_train=encoded_train.detach().cpu().numpy()
import pandas as pd
encoded_train=pd.DataFrame(encoded_train)
encoded_train["AACombo"]=labels
encoded_train.to_csv("UBE_VAE/UBE_homologues_representationz30.csv")


# embed new features?
### load test data ###
testdata_filename="UBE_VAE/UBE_mutants.fasta"
test_ids = []
test_labels = []
test_seqs = []
for record in SeqIO.parse(testdata_filename, "fasta"):
        test_ids.append(record.id)
        test_labels.append(record.description)
        test_seqs.append(np.array([aa1_to_index[aa] for aa in str(record.seq).upper().replace('.', '-')]))

test_seqs = torch.from_numpy(np.vstack(test_seqs))
test_labels = np.array(test_labels)
test_embedded=model.embedding(test_seqs) 

#### make to numpy array
test_embedded.cpu().detach().numpy()
np.save("UBE_VAE/UBE_variants_AACombo_VAEz30_embedded.npy", test_embedded.cpu().detach().numpy())
### make inx_to com and Com_to_ind
ind_to_combo={k:v for k, v in zip(range(len(test_labels)), test_labels)}
combo_to_ind={v:k for k, v in zip(range(len(test_labels)), test_labels)}
import pickle
with open('UBE_VAE/UBE_variants_AACombo_VAEz30_IndexToCombo.pkl', 'wb') as handle:
    pickle.dump(ind_to_combo, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('UBE_VAE/UBE_variants_AACombo_VAEz30_ComboToIndex.pkl', 'wb') as handle:
    pickle.dump(combo_to_ind, handle, protocol=pickle.HIGHEST_PROTOCOL)
