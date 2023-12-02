### ls /home/athchu/.cache/torch/hub/checkpoints/
### cd /data/users/athchu/efficient-evolution
### conda activate /data/users/athchu/.conda/envs/efficient-evolution
import torch
import esm
from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer
import sys
sys.path.append("path_to/efficient-evolution")
import pandas as pd
from Bio import SeqIO

### get prob from models ###
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.cuda(1)
batch_converter = alphabet.get_batch_converter()

def label_row(row, sequence, token_probs, alphabet, offset_idx):
    mut_scores=0
    wt_scores=0
    muts=row.split("_")
    for m in muts:
        wt, idx, mt = m[0], int(m[1:-1]) - offset_idx, m[-1]
        wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)
        mut_scores+=token_probs[0, 1 + idx, mt_encoded]
        wt_scores+=token_probs[0, 1 + idx, wt_encoded]
    return mut_scores, wt_scores

def generate_scores(aa, df, token_probs, alphabet):
	mut_scores=[]
	WT_scores=[]
	for i in df.index:
		row=df.loc[i]
		row=row["variants"]
		if row !="WT":
			mut_S, wt_S=label_row(row, aa, token_probs, alphabet, 0)
			mut_scores.append(mut_S.item())
			WT_scores.append(wt_S.item())
		else:
			mut_S, wt_S=0, 0 
			mut_scores.append(mut_S)
			WT_scores.append(wt_S)
	return mut_scores, WT_scores

### SaCas91296
aa=SeqIO.read("input_data/SaCas91296_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/SaCas91296_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm2_t33_650M_UR50D_mut"]=mut_scores
df["esm2_t33_650M_UR50D_WT"]=WT_scores
df.to_csv("SaCas91296_esm2_t33_650M_UR50D_WT_variants.csv")

### SaCAs98000
aa=SeqIO.read("input_data/SaCas98000_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/SaCas98000_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm2_t33_650M_UR50D_mut"]=mut_scores
df["esm2_t33_650M_UR50D_WT"]=WT_scores
df.to_csv("SaCas9800_esm2_t33_650M_UR50D_WT_variants.csv")

### SpCas9 
aa=SeqIO.read("input_data/SpWT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/Sp_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm2_t33_650M_UR50D_mut"]=mut_scores
df["esm2_t33_650M_UR50D_WT"]=WT_scores
df.to_csv("Sp_esm2_t33_650M_UR50D_WT_variants.csv")

### evoAPOBEC ###
aa=SeqIO.read("input_data/evoAPOBEC_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/evoAPOBEC_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm2_t33_650M_UR50D_mut"]=mut_scores
df["esm2_t33_650M_UR50D_WT"]=WT_scores
df.to_csv("evoAPOBEC_esm2_t33_650M_UR50D_WT_variants.csv")

### PhoQ
aa=SeqIO.read("input_data/PhoQ_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/PhoQ_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm2_t33_650M_UR50D_mut"]=mut_scores
df["esm2_t33_650M_UR50D_WT"]=WT_scores
df.to_csv("PhoQ_esm2_t33_650M_UR50D_WT_variants.csv")

### GFP
aa=SeqIO.read("input_data/GFP_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/GFP_variants.csv'
df = pd.read_csv(fname, delimiter=',')
df=df.rename(columns={"AACombo": "variants"})
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm2_t33_650M_UR50D_mut"]=mut_scores
df["esm2_t33_650M_UR50D_WT"]=WT_scores
df.to_csv("GFP_esm2_t33_650M_UR50D_WT_variants.csv")

### PABP
aa=SeqIO.read("input_data/PABP_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/PABP_variants.csv'
df = pd.read_csv(fname, delimiter=',')
df=df.rename(columns={"AACombo": "variants"})
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm2_t33_650M_UR50D_mut"]=mut_scores
df["esm2_t33_650M_UR50D_WT"]=WT_scores
df.to_csv("PABP_esm2_t33_650M_UR50D_WT_variants.csv")

### UBE
aa=SeqIO.read("input_data/UBE_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/UBE_variants.csv'
df = pd.read_csv(fname, delimiter=',')
df=df.rename(columns={"AACombo": "variants"})
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm2_t33_650M_UR50D_mut"]=mut_scores
df["esm2_t33_650M_UR50D_WT"]=WT_scores
df.to_csv("UBE_esm2_t33_650M_UR50D_WT_variants.csv")



### do it again with ESM-1v ###
model, alphabet = esm.pretrained.esm1v_t33_650M_UR90S_5()
model = model.cuda(1)
batch_converter = alphabet.get_batch_converter()

### SaCas91296
aa=SeqIO.read("input_data/SaCas91296_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/SaCas91296_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1v_t33_650M_UR90S_5_mut"]=mut_scores
df["esm1v_t33_650M_UR90S_5_WT"]=WT_scores
df.to_csv("SaCas91296_esm1v_t33_650M_UR90S_5_WT_variants.csv")

### SaCAs98000
aa=SeqIO.read("input_data/SaCas98000_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/SaCas98000_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1v_t33_650M_UR90S_5_mut"]=mut_scores
df["esm1v_t33_650M_UR90S_5_WT"]=WT_scores
df.to_csv("SaCas9800_esm1v_t33_650M_UR90S_5_WT_variants.csv")

### SpCas9 
aa=SeqIO.read("input_data/SpWT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/Sp_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1v_t33_650M_UR90S_5_mut"]=mut_scores
df["esm1v_t33_650M_UR90S_5_WT"]=WT_scores
df.to_csv("Sp_esm1v_t33_650M_UR90S_5_WT_variants.csv")

### evoAPOBEC ###
aa=SeqIO.read("input_data/evoAPOBEC_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/evoAPOBEC_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1v_t33_650M_UR90S_5_mut"]=mut_scores
df["esm1v_t33_650M_UR90S_5_WT"]=WT_scores
df.to_csv("evoAPOBEC_esm1v_t33_650M_UR90S_5_WT_variants.csv")

### PhoQ
aa=SeqIO.read("input_data/PhoQ_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/PhoQ_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1v_t33_650M_UR90S_5_mut"]=mut_scores
df["esm1v_t33_650M_UR90S_5_WT"]=WT_scores
df.to_csv("PhoQ_esm1v_t33_650M_UR90S_5_WT_variants.csv")

### GFP
aa=SeqIO.read("input_data/GFP_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/GFP_variants.csv'
df = pd.read_csv(fname, delimiter=',')
df=df.rename(columns={"AACombo": "variants"})
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1v_t33_650M_UR90S_5_mut"]=mut_scores
df["esm1v_t33_650M_UR90S_5_WT"]=WT_scores
df.to_csv("GFP_esm1v_t33_650M_UR90S_5_WT_variants.csv")

### PABP
aa=SeqIO.read("input_data/PABP_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/PABP_variants.csv'
df = pd.read_csv(fname, delimiter=',')
df=df.rename(columns={"AACombo": "variants"})
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1v_t33_650M_UR90S_5_mut"]=mut_scores
df["esm1v_t33_650M_UR90S_5_WT"]=WT_scores
df.to_csv("PABP_esm1v_t33_650M_UR90S_5_WT_variants.csv")

### UBE
aa=SeqIO.read("input_data/UBE_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/UBE_variants.csv'
df = pd.read_csv(fname, delimiter=',')
df=df.rename(columns={"AACombo": "variants"})
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1v_t33_650M_UR90S_5_mut"]=mut_scores
df["esm1v_t33_650M_UR90S_5_WT"]=WT_scores
df.to_csv("UBE_esm1v_t33_650M_UR90S_5_WT_variants.csv")

### do it again with ESM-1b ###
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
model = model.cuda(1)
batch_converter = alphabet.get_batch_converter()

### SaCas91296
aa=SeqIO.read("input_data/SaCas91296_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/SaCas91296_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1b_t33_650M_UR50S_mut"]=mut_scores
df["esm1b_t33_650M_UR50S_WT"]=WT_scores
df.to_csv("SaCas91296_esm1b_t33_650M_UR50S_WT_variants.csv")

### SaCAs98000
aa=SeqIO.read("input_data/SaCas98000_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/SaCas98000_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1b_t33_650M_UR50S_mut"]=mut_scores
df["esm1b_t33_650M_UR50S_WT"]=WT_scores
df.to_csv("SaCas9800_esm1b_t33_650M_UR50S_WT_variants.csv")

### SpCas9 
aa=SeqIO.read("input_data/SpWT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/Sp_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1b_t33_650M_UR50S_mut"]=mut_scores
df["esm1b_t33_650M_UR50S_WT"]=WT_scores
df.to_csv("Sp_esm1b_t33_650M_UR50S_WT_variants.csv")

### evoAPOBEC ###
aa=SeqIO.read("input_data/evoAPOBEC_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/evoAPOBEC_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1b_t33_650M_UR50S_mut"]=mut_scores
df["esm1b_t33_650M_UR50S_WT"]=WT_scores
df.to_csv("evoAPOBEC_esm1b_t33_650M_UR50S_WT_variants.csv")

### PhoQ
aa=SeqIO.read("input_data/PhoQ_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/PhoQ_variants.csv'
df = pd.read_csv(fname, delimiter=',')
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1b_t33_650M_UR50S_mut"]=mut_scores
df["esm1b_t33_650M_UR50S_WT"]=WT_scores
df.to_csv("PhoQ_esm1b_t33_650M_UR50S_WT_variants.csv")

### GFP
aa=SeqIO.read("input_data/GFP_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/GFP_variants.csv'
df = pd.read_csv(fname, delimiter=',')
df=df.rename(columns={"AACombo": "variants"})
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1b_t33_650M_UR50S_mut"]=mut_scores
df["esm1b_t33_650M_UR50S_WT"]=WT_scores
df.to_csv("GFP_esm1b_t33_650M_UR50S_WT_variants.csv")

### PABP
aa=SeqIO.read("input_data/PABP_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/PABP_variants.csv'
df = pd.read_csv(fname, delimiter=',')
df=df.rename(columns={"AACombo": "variants"})
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1b_t33_650M_UR50S_mut"]=mut_scores
df["esm1b_t33_650M_UR50S_WT"]=WT_scores
df.to_csv("PABP_esm1b_t33_650M_UR50S_WT_variants.csv")

### UBE
aa=SeqIO.read("input_data/UBE_WT.fasta", "fasta")
aa=str(aa.seq)
data = [("protein1", aa),]
batch_labels, batch_strs, batch_tokens = batch_converter(data)
with torch.no_grad():
    token_probs = torch.log_softmax(model(batch_tokens.cuda(1))["logits"], dim=-1)

fname='input_data/UBE_variants.csv'
df = pd.read_csv(fname, delimiter=',')
df=df.rename(columns={"AACombo": "variants"})
mut_scores, WT_scores=generate_scores(aa, df, token_probs, alphabet)
df["esm1b_t33_650M_UR50S_mut"]=mut_scores
df["esm1b_t33_650M_UR50S_WT"]=WT_scores
df.to_csv("UBE_esm1b_t33_650M_UR50S_WT_variants.csv")

