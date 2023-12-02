import numpy as np
import os
import sys
import pandas as pd
from sklearn.cluster import KMeans
import warnings
import pickle
import copy
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
import matplotlib.pyplot as plt
import seaborn as sn

### cluster_sample ###
def cluster_extract_sampling_prob(data_zeroshot, n_clusters, ComboToIndex, IndextoCombo, AACombo, features, zeroshot="predictor"):
    num_first_round=96
    if 'cluster_id' not in data_zeroshot.columns:
        ### Run k_mean
        kmeans = KMeans(n_clusters=n_clusters).fit(features)
        ### run Clustering
        cluster_labels = kmeans.labels_
        ### map features to AACombo
        clusters=pd.DataFrame(list(zip(AACombo, kmeans.labels_)), columns=["AACombo", "cluster_id"])
        ### 
        data_zeroshot = pd.merge(data_zeroshot, clusters, how="left", on="AACombo")
    ### store selected samples with sequential order
    Cluster_list=[]
    #  store selected samples according to the cluster they belong to
    # first round sampling
    ### sample based on probability of cluster # here ==1 atm ###
    ### rewrite with zeroshot
    Fit = [[] for _ in range(n_clusters)]
    SEQ = [[] for _ in range(n_clusters)]
    SEQ_index = [[] for _ in range(n_clusters)]
    num = 0
    Prob = np.ones([n_clusters]) / n_clusters
    Cluster_list=[0]*25
    while 0 in Cluster_list:
        Cluster_list=list(np.random.choice(np.arange(0, n_clusters), size=num_first_round, replace=True, p=Prob))
        Cluster_list=[Cluster_list.count(i) for i in range(0, n_clusters)]
    for i in range(0, n_clusters):
        num_var=Cluster_list[i]
        if num_var <1:
            continue
        else:
            temp_tab=data_zeroshot.loc[data_zeroshot["cluster_id"]==i]
            if temp_tab.shape[0]>num_var:
                temp_tab=temp_tab.sample(num_var, axis=0) ### only works with we are sure num_var << num. of rows
            elif temp_tab.shape[0]==0:
                break
            SEQ[i].extend(temp_tab["AACombo"])
            Fit[i].extend(temp_tab[zeroshot])
            SEQ_index[i].extend([ComboToIndex[c] for c in SEQ[i]])
    #### update Prob so that higher fitness cluster has a higher drawing prob.
    Mean_Fit = np.asarray([np.nanmean(np.asarray(Fit[i])) for i in range(n_clusters)])
    min_val=np.nanmin(Mean_Fit)
    Mean_Fit[np.where(np.isnan(Mean_Fit))[0]] = min_val
    # min_max normalise? 
    Mean_Fit =(Mean_Fit - min(Mean_Fit))/(max(Mean_Fit) - min(Mean_Fit))
    ### correct for neg values
    #### now what to do with Mean_Fit that seek less negative values? ###
    Prob = Mean_Fit / np.sum(Mean_Fit)
    if len(Prob[np.where(Prob<=0)])>0:
        Prob[np.where(Prob<=0)] = np.min(Prob[np.where(Prob>0)])
        Prob=Prob/sum(Prob)
    #### to correct for some cluster has too small probability of getting samples
    if max(Prob)/min(Prob) >10:
        Prob[np.where(Prob==min(Prob))] = max(Prob)/10
        Prob=Prob/sum(Prob)
    return Prob, data_zeroshot

def compute_uncertainty(Prob, data_zeroshot, features, size, zeroshot="predictor"):
    """sampling_subcluster_priority """
    seed=100
    #acquisition="UCB"
    sampling_para=4
    ### fit guassian Kde
    #if acquisition in ['UCB', 'epsilon','Thompson']:
    pred_mean_list=[]
    pred_std_list=[]
    kernel = 1 * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e2))
    ### cannot use all zeroshot data
    ### use 5% of data? and run 4 batch?
    for i in range(4): 
        temp_tab=data_zeroshot.sample(n=100)
        temp_tab=temp_tab.dropna(subset=["cluster_id", zeroshot])
        temp_tab["feature_index"]=temp_tab["feature_index"].astype(int)
        X_GP=features[temp_tab["feature_index"]]
        Y_GP=temp_tab[zeroshot]
        Y_GP=np.asarray(Y_GP)
        regr = GaussianProcessRegressor(random_state=seed, kernel=kernel, n_restarts_optimizer=9)
        regr.fit(X_GP, Y_GP)
        pred_mean, pred_std = regr.predict(features, return_std=True)
        pred_mean_list.append(pred_mean)
        pred_std_list.append(pred_std)
    pred_mean=np.asarray(pred_mean_list)
    pred_mean=np.mean(pred_mean, axis=0)
    pred_std=np.asarray(pred_std_list)
    pred_std=np.mean(pred_std, axis=0)
    pred_total=pred_mean + pred_std * np.sqrt(sampling_para)
    pred=pd.DataFrame(list(zip(AACombo, pred_mean, pred_std, pred_total)), columns=["AACombo", "pred_mean", "pred_std", "pred_total"])
    data_zeroshot = pd.merge(data_zeroshot, pred, how="outer", on="AACombo")
    return data_zeroshot


### draw sample based on updated PROB
def extract_CLADE_aa(data_zeroshot, n_clusters, Prob, size, zeroshot, orders=False):
    Cluster_list=[0]*25
    n=0
    while 0 in Cluster_list and n<25:
        """iterate certain time if it is impossible to sample certain clusters """
        Cluster_list=list(np.random.choice(np.arange(0, n_clusters), size=size, replace=True, p=Prob))
        Cluster_list=[Cluster_list.count(i) for i in range(0, n_clusters)]
        n+=1
### take the top (least negative variants) for each clusters ###
### a higher (less negative) zero-shot score means greater confidence that that variant will be functional. ###
    select_AACombo=[]
    for i in range(0, n_clusters):
        num_var=Cluster_list[i]
        temp_tab=data_zeroshot.loc[data_zeroshot["cluster_id"]==i]
        temp_tab=temp_tab.sort_values(by=zeroshot, ascending=orders)
        temp_tab=temp_tab[0:num_var]
        select_AACombo.extend(list(temp_tab["AACombo"]))
    return select_AACombo


### MSA 
### load in data
zero_shot="/media/achu/新增磁碟區/MLDE/2021_11_09_new_MLDE/SpCas9_586-1255/SpCas9_R76_Q110_K263_E338_T339_Q341_K418_MSA_ZeroShotPreds.csv"
zeroshot="esm_msa1_t12_100M_UR50S-Naive-ColumnMasked"
data_zeroshot = pd.read_csv(zero_shot)
data_zeroshot = data_zeroshot.rename({"Combo":"AACombo"}, axis=1)
#### MSA zeor-shot has 1060R missing 
data_zeroshot["AACombo"]= data_zeroshot["AACombo"] + "R"

#### import encoding
encoding_file="data/clustered_sampling/SpCas9_georgiev_UnNormalized.npy"
ComboToIndex_file="SpCas9_georgiev_ComboToIndex.pkl"
IndextoCombo_file="SpCas9_georgiev_IndexToCombo.pkl"
features=np.load(encoding_file)
ComboToIndex=pickle.load(open(ComboToIndex_file, "rb"))
IndextoCombo=pickle.load(open(IndextoCombo_file, "rb"))

### add index info to table ###
data_zeroshot =data_zeroshot.loc[data_zeroshot["AACombo"].isin(list(ComboToIndex.keys()))]
data_zeroshot = pd.merge(data_zeroshot,  
pd.DataFrame(zip(ComboToIndex.keys(), ComboToIndex.values()), columns=["AACombo", "feature_index"]),
how="right", on="AACombo")
data_zeroshot["feature_index"]=data_zeroshot["feature_index"].astype(int)
data_zeroshot= data_zeroshot.drop_duplicates()
data_zeroshot= data_zeroshot.dropna(subset=["feature_index"])

### normalise zero-shot predictor for Guassian regressor, assume all are negative values 
data_zeroshot["norm"] = (data_zeroshot[zeroshot] - np.nanmin(data_zeroshot[zeroshot]))/(np.nanmax(data_zeroshot[zeroshot])- np.nanmin(data_zeroshot[zeroshot]))
data_zeroshot["predictor"] = (data_zeroshot[zeroshot] - np.mean(data_zeroshot[zeroshot])) /np.std(data_zeroshot[zeroshot])
### make features
AACombo=[IndextoCombo[i] for i in range(0, len(IndextoCombo.keys()))]
if len(features.shape) == 3:
    features = np.reshape(features, [features.shape[0], features.shape[1] * features.shape[2]])

### make k-meansdistortions = []
import matplotlib.pyplot as plt

distortions = []
K = range(1,15)
for k in K:
    kmeanModel = KMeans(n_clusters=k)
    kmeanModel.fit(features)
    distortions.append(kmeanModel.inertia_)

plt.figure(figsize=(16,8))
plt.plot(K, distortions, 'bx-')
plt.axvline(x=10, color='k', linestyle='--')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('The Elbow Method showing the optimal k')
plt.show()
### keep all features, even without any zero shot prediction? ###
#features = features[data_zeroshot["feature_index"]]
Prob, data_zeroshot=cluster_extract_sampling_prob(data_zeroshot, 10, ComboToIndex, IndextoCombo, AACombo, features, "norm")

wkdir="data/clustered_sampling/SpCas9_MLDE_Fitness_input/"
outprefix="SpCas9_R76Q110K263E338T339Q341K418_CLADE_MSA_"
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 96, zeroshot, orders=False)
o=open(wkdir+outprefix+"96.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 48, zeroshot, orders=False)
o=open(wkdir+outprefix+"48.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 24, zeroshot, orders=False)
o=open(wkdir+outprefix+"24.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 12, zeroshot, orders=False)
o=open(wkdir+outprefix+"12.txt", "w")
o.write("\n".join(select_AACombo))
o.close()


#### add uncertainty
Prob, data_zeroshot=cluster_extract_sampling_prob(data_zeroshot, 10, ComboToIndex, IndextoCombo, AACombo, features, "predictor")
data_zeroshot=compute_uncertainty(Prob, data_zeroshot, features, 100, "predictor")
#### check sklearn results
import matplotlib.pyplot as plt
import seaborn as sn
sn.scatterplot(data=data_zeroshot, x="predictor", y="pred_mean")
plt.show()
sn.scatterplot(data=data_zeroshot, x="predictor", y=zeroshot)
plt.show()
sn.scatterplot(data=data_zeroshot, x=zeroshot, y="pred_total")
plt.show()
plt.hist(data_zeroshot["pred_std"], bins=50)
plt.show()
data_zeroshot.to_csv(wkdir+"SpCas9_R76Q110K263E338T339Q341K418_MSA_clustered_ZeroShotPreds.csv")

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 96, "pred_total", orders=False)
o=open(wkdir+outprefix+"96_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 48, "pred_total", orders=False)
o=open(wkdir+outprefix+"48_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 24, "pred_total", orders=False)
o=open(wkdir+outprefix+"24_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 12, "pred_total", orders=False)
o=open(wkdir+outprefix+"12_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

### EVmutation
### load in data
zero_shot="data/clustered_sampling/SPCAS9_CombiSEAL_R661_Q695_K848_E923_T924_Q926_K1003.csv"
zeroshot="EvMutation"
data_zeroshot = pd.read_csv(zero_shot)
data_zeroshot = data_zeroshot.rename({"Combo":"AACombo"}, axis=1)
#### MSA zeor-shot has 1060R missing 
data_zeroshot["AACombo"]= data_zeroshot["AACombo"] + "R"
data_zeroshot= data_zeroshot.dropna(subset=["feature_index"])

#### import encoding
encoding_file="SpCas9_georgiev_UnNormalized.npy"
ComboToIndex_file="SpCas9_georgiev_ComboToIndex.pkl"
IndextoCombo_file="SpCas9_georgiev_IndexToCombo.pkl"
features=np.load(encoding_file)
ComboToIndex=pickle.load(open(ComboToIndex_file, "rb"))
IndextoCombo=pickle.load(open(IndextoCombo_file, "rb"))

### add index info to table ###
data_zeroshot =data_zeroshot.loc[data_zeroshot["AACombo"].isin(list(ComboToIndex.keys()))]
data_zeroshot = pd.merge(data_zeroshot,  
pd.DataFrame(zip(ComboToIndex.keys(), ComboToIndex.values()), columns=["AACombo", "feature_index"]),
how="right", on="AACombo")
data_zeroshot["feature_index"]=data_zeroshot["feature_index"].astype(int)
data_zeroshot= data_zeroshot.drop_duplicates()
### keep entry that are in the features for clustering even it doesn't has the zeroshot value
data_zeroshot= data_zeroshot.dropna(subset=["feature_index"])

### normalise zero-shot predictor for Guassian regressor, assume all are negative values 
data_zeroshot["norm"] = (data_zeroshot[zeroshot] - np.nanmin(data_zeroshot[zeroshot]))/(np.nanmax(data_zeroshot[zeroshot])- np.nanmin(data_zeroshot[zeroshot]))
data_zeroshot["predictor"] = (data_zeroshot[zeroshot] - np.mean(data_zeroshot[zeroshot])) /np.std(data_zeroshot[zeroshot])
### make features
AACombo=[IndextoCombo[i] for i in range(0, len(IndextoCombo.keys()))]
if len(features.shape) == 3:
    features = np.reshape(features, [features.shape[0], features.shape[1] * features.shape[2]])

Prob, data_zeroshot=cluster_extract_sampling_prob(data_zeroshot, 10, ComboToIndex, IndextoCombo, AACombo, features, "norm")

wkdir="data/clustered_sampling/SpCas9_MLDE_Fitness_input/"
outprefix="SpCas9_R661Q695K848E923T924Q926K1003_CLADE_Evmutation_"
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 96, zeroshot, orders=False)
o=open(wkdir+outprefix+"96.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 48, zeroshot, orders=False)
o=open(wkdir+outprefix+"48.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 24, zeroshot, orders=False)
o=open(wkdir+outprefix+"24.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 12, zeroshot, orders=False)
o=open(wkdir+outprefix+"12.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

#### add uncertainty
Prob, data_zeroshot=cluster_extract_sampling_prob(data_zeroshot, 10, ComboToIndex, IndextoCombo, AACombo, features, "predictor")
data_zeroshot=compute_uncertainty(Prob, data_zeroshot, features, 100, "predictor")
#### check sklearn results
import matplotlib.pyplot as plt
import seaborn as sn
sn.scatterplot(data=data_zeroshot, x="predictor", y="pred_mean")
plt.show()
sn.scatterplot(data=data_zeroshot, x="predictor", y=zeroshot)
plt.show()
plt.hist(data_zeroshot["pred_std"], bins=50)
plt.show()

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 96, "pred_total", orders=False)
o=open(wkdir+outprefix+"96_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 48, "pred_total", orders=False)
o=open(wkdir+outprefix+"48_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 24, "pred_total", orders=False)
o=open(wkdir+outprefix+"24_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 12, "pred_total", orders=False)
o=open(wkdir+outprefix+"12_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

data_zeroshot.to_csv(wkdir + "SpCas9_R661Q695K848E923T924Q926K1003_EvMutation_clustered_ZeroShotPreds.csv")

### EvoEF
### reconstructu EvoEF full table
zero_shot="data/clustered_sampling/SpCas9_R661Q695K848E923T924Q926K1003R1060_EvoEF.csv"
zeroshot="DDG_stab_1"
data_zeroshot = pd.read_csv(zero_shot)

#### import encoding
encoding_file="SpCas9_georgiev_UnNormalized.npy"
ComboToIndex_file="SpCas9_georgiev_ComboToIndex.pkl"
IndextoCombo_file="SpCas9_georgiev_IndexToCombo.pkl"
features=np.load(encoding_file)
ComboToIndex=pickle.load(open(ComboToIndex_file, "rb"))
IndextoCombo=pickle.load(open(IndextoCombo_file, "rb"))

### add index info to table ###
data_zeroshot = pd.merge(data_zeroshot,  
pd.DataFrame(zip(ComboToIndex.keys(), ComboToIndex.values()), columns=["AACombo", "feature_index"]),
how="right", on="AACombo")

### DDG _stab is a bit different, the more negative it is the better it is 
### so * -1 to reverse the scale
data_zeroshot["DDG_stab_1"]=data_zeroshot["DDG_stab"]  * -1
sn.scatterplot(data=data_zeroshot, x="DDG_stab", y="DDG_stab_1")
plt.show()

### normalise zero-shot predictor for Guassian regressor, assume all are negative values 
data_zeroshot["norm"] = (data_zeroshot[zeroshot] - np.nanmin(data_zeroshot[zeroshot]))/(np.nanmax(data_zeroshot[zeroshot])- np.nanmin(data_zeroshot[zeroshot]))
data_zeroshot["predictor"] = (data_zeroshot[zeroshot] - np.mean(data_zeroshot[zeroshot])) /np.std(data_zeroshot[zeroshot])
### make features
AACombo=[IndextoCombo[i] for i in range(0, len(IndextoCombo.keys()))]
if len(features.shape) == 3:
    features = np.reshape(features, [features.shape[0], features.shape[1] * features.shape[2]])


Prob, data_zeroshot=cluster_extract_sampling_prob(data_zeroshot, 10, ComboToIndex, IndextoCombo, AACombo, features, "norm")
wkdir="data/clustered_sampling/SpCas9_MLDE_Fitness_input/"
outprefix="SpCas9_R661Q695K848E923T924Q926K1003R1060_CLADE_EvoEF_"

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 96, zeroshot, orders=False)
o=open(wkdir+outprefix+"96.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 48, zeroshot, orders=False)
o=open(wkdir+outprefix+"48.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 24, zeroshot, orders=False)
o=open(wkdir+outprefix+"24.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 12, zeroshot, orders=False)
o=open(wkdir+outprefix+"12.txt", "w")
o.write("\n".join(select_AACombo))
o.close()


#### add uncertainty
Prob, data_zeroshot=cluster_extract_sampling_prob(data_zeroshot, 10, ComboToIndex, IndextoCombo, AACombo, features, "predictor")
data_zeroshot=compute_uncertainty(Prob, data_zeroshot, features, 100, "predictor")
#### check sklearn results
import matplotlib.pyplot as plt
import seaborn as sn
sn.scatterplot(data=data_zeroshot, x="predictor", y="pred_mean")
plt.show()
sn.scatterplot(data=data_zeroshot, x="predictor", y=zeroshot)
plt.show()
plt.hist(data_zeroshot["pred_std"], bins=50)
plt.show()

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 96, "pred_total", orders=False)
o=open(wkdir+outprefix+"96_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 48, "pred_total", orders=False)
o=open(wkdir+outprefix+"48_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 24, "pred_total", orders=False)
o=open(wkdir+outprefix+"24_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 12, "pred_total", orders=False)
o=open(wkdir+outprefix+"12_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

data_zeroshot.to_csv(wkdir + "SpCas9_R661Q695K848E923T924Q926K1003R1060_CLADE_EvoEF_clustered_ZeroShotPreds.csv")

### arDCA
### reconstructu EvoEF full table
zero_shot="data/clustered_sampling/SpCas9_R661Q695K848E923T924Q926K1003R1060_arDCA.csv"
zeroshot="arDCA_1"
data_zeroshot = pd.read_csv(zero_shot)

#### import encoding
encoding_file="SpCas9_georgiev_UnNormalized.npy"
ComboToIndex_file="SpCas9_georgiev_ComboToIndex.pkl"
IndextoCombo_file="SpCas9_georgiev_IndexToCombo.pkl"
features=np.load(encoding_file)
ComboToIndex=pickle.load(open(ComboToIndex_file, "rb"))
IndextoCombo=pickle.load(open(IndextoCombo_file, "rb"))

### add index info to table ###
data_zeroshot = pd.merge(data_zeroshot,  
pd.DataFrame(zip(ComboToIndex.keys(), ComboToIndex.values()), columns=["AACombo", "feature_index"]),
how="right", on="AACombo")

### DDG _stab is a bit different, the more negative it is the better it is 
### so * -1 to reverse the scale
data_zeroshot["arDCA_1"]=data_zeroshot["arDCA"]  * -1
sn.scatterplot(data=data_zeroshot, x="arDCA", y="arDCA_1")
plt.show()

### normalise zero-shot predictor for Guassian regressor, assume all are negative values 
data_zeroshot["norm"] = (data_zeroshot[zeroshot] - np.nanmin(data_zeroshot[zeroshot]))/(np.nanmax(data_zeroshot[zeroshot])- np.nanmin(data_zeroshot[zeroshot]))
data_zeroshot["predictor"] = (data_zeroshot[zeroshot] - np.mean(data_zeroshot[zeroshot])) /np.std(data_zeroshot[zeroshot])
### make features
AACombo=[IndextoCombo[i] for i in range(0, len(IndextoCombo.keys()))]
if len(features.shape) == 3:
    features = np.reshape(features, [features.shape[0], features.shape[1] * features.shape[2]])


Prob, data_zeroshot=cluster_extract_sampling_prob(data_zeroshot, 10, ComboToIndex, IndextoCombo, AACombo, features, "norm")
wkdir="data/clustered_sampling/SpCas9_MLDE_Fitness_input/"
outprefix="SpCas9_R661Q695K848E923T924Q926K1003R1060_CLADE_arDCA_"

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 96, zeroshot, orders=False)
o=open(wkdir+outprefix+"96.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 48, zeroshot, orders=False)
o=open(wkdir+outprefix+"48.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 24, zeroshot, orders=False)
o=open(wkdir+outprefix+"24.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 12, zeroshot, orders=False)
o=open(wkdir+outprefix+"12.txt", "w")
o.write("\n".join(select_AACombo))
o.close()


#### add uncertainty
Prob, data_zeroshot=cluster_extract_sampling_prob(data_zeroshot, 10, ComboToIndex, IndextoCombo, AACombo, features, "predictor")
data_zeroshot=compute_uncertainty(Prob, data_zeroshot, features, 100, "predictor")
#### check sklearn results
import matplotlib.pyplot as plt
import seaborn as sn
sn.scatterplot(data=data_zeroshot, x="predictor", y="pred_mean")
plt.show()
sn.scatterplot(data=data_zeroshot, x="predictor", y=zeroshot)
plt.show()
plt.hist(data_zeroshot["pred_std"], bins=50)
plt.show()

select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 96, "pred_total", orders=False)
o=open(wkdir+outprefix+"96_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 48, "pred_total", orders=False)
o=open(wkdir+outprefix+"48_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 24, "pred_total", orders=False)
o=open(wkdir+outprefix+"24_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()
select_AACombo=extract_CLADE_aa(data_zeroshot, 10, Prob, 12, "pred_total", orders=False)
o=open(wkdir+outprefix+"12_gs.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

data_zeroshot.to_csv(wkdir + "SpCas9_R661Q695K848E923T924Q926K1003R1060_CLADE_arDCA_clustered_ZeroShotPreds.csv")





### random 
select_AACombo=list(data_zeroshot.sample(96, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep1_random_96.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(48, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep1_random_48.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(24, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep1_random_24.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(12, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep1_random_12.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(96, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep2_random_96.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(48, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep2_random_48.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(24, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep2_random_24.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(12, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep2_random_12.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(96, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep3_random_96.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(48, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep3_random_48.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(24, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep3_random_24.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

select_AACombo=list(data_zeroshot.sample(12, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_rep3_random_12.txt", "w")
o.write("\n".join(select_AACombo))
o.close()

### fix 20% of data as test data ###
test_AACombo=list(data_zeroshot.sample(190, axis=0)["AACombo"])
o=open(wkdir+"SpCas9_R661Q695K848E923T924Q926K1003R1060_fixed_test_set.txt", "w")
o.write("\n".join(test_AACombo))
o.close()
### no need to remove test set becuase these are zero shot predictions ###

### make fitness.csv for MLDE
Fitness=pd.read_csv("data/MLDE_train_test_datasets/SpCas9_fitness_norm.csv")
indir="data/clustered_sampling/SpCas9_MLDE_Fitness_input/"
outdir="data/clustered_sampling/SpCas9_MLDE_Fitness_input/220705_SpCas9_R661Q695K848E923T924Q926K1003R1060_Fitnesscsv/"

infiles=os.listdir(indir)
infiles=[i for i in infiles if "SpCas9" in i]
infiles=[i for i in infiles if "txt" in i]
infiles=[i for i in infiles if "fixed" not in i]
infiles=[i for i in infiles if "96" in i or "48" in i or "24" in i or "12" in i]


for f in infiles:
    infile=f
    select_AACombo=pd.read_csv(indir+infile, header=None)
    select_AACombo=list(select_AACombo[0])
    temp_tab=Fitness[Fitness["AACombo"].isin(select_AACombo)]
    o=temp_tab[["AACombo", "Sg5"]]
    o=o.rename({"Sg5":"Fitness"}, axis=1)
    o=o.dropna(axis=0, how='any')
    outfile=infile.replace(".txt", "_Sg5.csv")
    o.to_csv(outdir+outfile, index=False)
    o=temp_tab[["AACombo", "Sg8"]]
    o=o.rename({"Sg8":"Fitness"}, axis=1)
    o=o.dropna(axis=0, how='any')
    outfile=infile.replace(".txt", "_Sg8.csv")
    o.to_csv(outdir+outfile, index=False)
