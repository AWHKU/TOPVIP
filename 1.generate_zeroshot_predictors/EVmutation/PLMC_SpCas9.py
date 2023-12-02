#### conda activate protein_fitness_prediction
####
import argparse
from multiprocessing import Process, JoinableQueue
from multiprocessing import set_start_method
import os
import pandas as pd
import sys
sys.path.append('path_to/combining-evolutionary-and-assay-labelled-data/src')

from predictors import base_predictors, ev_predictors, hmm_predictors, onehot_predictors
from predictors import unirep_predictors, esm_predictors, vae_predictors
from predictors.base_predictors import BoostingPredictor, JointPredictor
from utils.metric_utils import spearman, topk_mean, r2, hit_rate, aucroc, ndcg
from utils.io_utils import load_data_split, get_wt_log_fitness, get_log_fitness_cutoff
from utils.data_utils import dict2str

predictor_cls = [ev_predictors.EVPredictor, onehot_predictors.OnehotRidgePredictor]
predictor = JointPredictor(dataset_name, predictor_cls, predictor_name, **predictor_params)


from utils import parse_vars, merge_dfs
from evaluate_multiprocessing import run_from_queue
from predictors import get_predictor_names 
dataset_name="SPCAS9_586-1255"

#### make ev_predictor
from utils import seqs_to_onehot, get_wt_seq, read_fasta, seq2effect, mutant2seq
from predictors.base_predictors import BaseRegressionPredictor
from couplings_model import CouplingsModel

model_name="uniref100"
m_path="data/zeroshot_prediction/EVcoupling_alignments/SaCas9/J7RUA5_b0.1.model"
couplings_model = CouplingsModel(m_path)
### begin couplings_model.py 
precision="float32"
file_format="plmc_v2"
wtseqs, wtids = read_fasta("data/zeroshot_prediction/EVcoupling/SPCAS9_e40.fa", return_ids=True)

couplings_model.index_list
### so that all the valid columns in the alignment
if '/' in wtids[0]:
    couplings_model.offset = int(wtids[0].split('/')[-1].split('-')[0])
else:
    couplings_model.offset = 1


expected_wt = wtseqs[0]
for pf, pm in couplings_model.index_map.items():
            if expected_wt[pf-couplings_model.offset] != couplings_model.target_seq[pm]:
                print(f'WT and model target seq mismatch at {pf}')
### al mismatch because target-seq is "-" ###

data=pd.read_csv("data/zeroshot_prediction/arDCA/SpCas9_fitness.csv")
X=seq2effect(data.seq.values, couplings_model, couplings_model.offset, ignore_gaps=False)
data['plmc']=X 

data.to_csv("SpCas9_EV.csv")
