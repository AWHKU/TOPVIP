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

### merge tables of empirical, EV, arDCA, seqdesign, cluster_id ###
### GFP
encoding_file="GFP_georgiev_embedding.npy"
ComboToIndex_file="GFP_georgiev_ComboToIndex.pkl"
EV_file="GFP_fitness_EV_rename.csv"
arDCA_file="data/clustered_sampling/GFP_arDCA.csv"
arDCA_input="data/clustered_sampling/GFP_arDCA_input.csv"
seqdesign_file="data/clustered_sampling/GFP_100_mutants.csv"
EV=pd.read_csv(EV_file)
### log_fitness, 'mutations', plmc

### process arDCA to get correct mutations 
arDCA=pd.read_csv(arDCA_file)
arDCA=arDCA.set_index("AACombo", drop=False)
arDCA_input=pd.read_csv(arDCA_input)
arDCA_input=arDCA_input.set_index('converted_mutations', drop=False)
arDCA=arDCA.join(arDCA_input)
arDCA=arDCA[['log_fitness', 'arDCA', 'mutations']]
arDCA=arDCA.set_index("mutations", drop=False)

### 'arDCA', 'AACombo'
seqdesign=pd.read_csv(seqdesign_file)
seqdesign=seqdesign.set_index("name", drop=False)
seqdesign=seqdesign.rename(columns={'mean': 'seqdesign'})
### 'name', 'mean'
df=EV[['log_fitness', 'mutations', 'plmc']]
df = df.rename(columns={'log_fitness': 'Fitness', 'mutations': 'AACombo'})
df=df.set_index("AACombo", drop=False)
df=df.join(arDCA[["arDCA",'log_fitness']])
df=df.join(seqdesign[["name", "seqdesign"]])

features=np.load(encoding_file)
ComboToIndex=pickle.load(open(ComboToIndex_file, "rb"))
AACombo=list(ComboToIndex.keys())
if len(features.shape) == 3:
    features = np.reshape(features, [features.shape[0], features.shape[1] * features.shape[2]])

distortions = []
K = range(1,30, 2)
for k in K:
    kmeanModel = KMeans(n_clusters=k)
    kmeanModel.fit(features)
    distortions.append(kmeanModel.inertia_)

plt.figure(figsize=(16,8))
plt.plot(K, distortions, 'bx-')
#plt.axvline(x=10, color='k', linestyle='--')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('The Elbow Method showing the optimal k')
plt.show()
### set cluster to 30
kmeans = KMeans(n_clusters=30).fit(features)
cluster_labels = kmeans.labels_
clusters=pd.DataFrame(list(zip(AACombo, kmeans.labels_)), columns=["AACombo", "cluster_id"])
clusters=clusters.set_index("AACombo", drop=False)
df=df.join(clusters["cluster_id"])
df["Fitness_n"]=(df["Fitness"].values - min(df["Fitness"].values))/ (max(df["Fitness"].values)-min(df["Fitness"].values))
df[['AACombo', 'Fitness', 'Fitness_n', 'plmc', 'arDCA','seqdesign', 'cluster_id']].\
to_csv("data/clustered_sampling/GFP_zeroshot.csv")

### PABP 
encoding_file="PABP_georgiev_embedding.npy"
ComboToIndex_file="PABP_georgiev_ComboToIndex.pkl"
EV_file="data/clustered_sampling/PABP_fitness_EV_rename.csv"
arDCA_file="data/clustered_sampling/PABP_arDCA.csv"
arDCA_input="data/clustered_sampling/PABP_arDCA_input.csv"
seqdesign_file="data/clustered_sampling/PABP_100_mutants.csv"
EV=pd.read_csv(EV_file)

### process arDCA to get correct mutations 
arDCA=pd.read_csv(arDCA_file)
arDCA=arDCA.set_index("AACombo", drop=False)
arDCA_input=pd.read_csv(arDCA_input)
arDCA_input=arDCA_input.set_index('converted_mutations', drop=False)
arDCA=arDCA.join(arDCA_input)
arDCA=arDCA[['log_fitness', 'arDCA', 'mutations']]
arDCA=arDCA.set_index("mutations", drop=False)

### 'arDCA', 'AACombo'
seqdesign=pd.read_csv(seqdesign_file)
seqdesign=seqdesign.set_index("name", drop=False)
seqdesign=seqdesign.rename(columns={'mean': 'seqdesign'})
### 'name', 'mean'
df=EV[['log_fitness', 'mutations', 'plmc']]
df = df.rename(columns={'log_fitness': 'Fitness', 'mutations': 'AACombo'})
df=df.set_index("AACombo", drop=False)
df=df.join(arDCA[["arDCA",'log_fitness']])
df=df.join(seqdesign[["name", "seqdesign"]])

features=np.load(encoding_file)
ComboToIndex=pickle.load(open(ComboToIndex_file, "rb"))
AACombo=list(ComboToIndex.keys())
if len(features.shape) == 3:
    features = np.reshape(features, [features.shape[0], features.shape[1] * features.shape[2]])

distortions = []
K = range(1,30, 2)
for k in K:
    kmeanModel = KMeans(n_clusters=k)
    kmeanModel.fit(features)
    distortions.append(kmeanModel.inertia_)

plt.figure(figsize=(16,8))
plt.plot(K, distortions, 'bx-')
#plt.axvline(x=10, color='k', linestyle='--')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('The Elbow Method showing the optimal k')
plt.show()
### set cluster to 30
kmeans = KMeans(n_clusters=27).fit(features)
cluster_labels = kmeans.labels_
clusters=pd.DataFrame(list(zip(AACombo, kmeans.labels_)), columns=["AACombo", "cluster_id"])
clusters=clusters.set_index("AACombo", drop=False)
df=df.join(clusters["cluster_id"])
df["Fitness_n"]=(df["Fitness"].values - min(df["Fitness"].values))/ (max(df["Fitness"].values)-min(df["Fitness"].values))
df[['AACombo', 'Fitness', 'Fitness_n', 'plmc', 'arDCA','seqdesign', 'cluster_id']].\
to_csv("data/clustered_sampling/PABP_zeroshot.csv")

### UBE ###
encoding_file="UBE_georgiev_embedding.npy"
ComboToIndex_file="UBE_georgiev_ComboToIndex.pkl"
EV_file="data/clustered_sampling/UBE_fitness_EV_rename.csv"
arDCA_file="data/clustered_sampling/UBE_arDCA.csv"
arDCA_input="data/clustered_sampling/UBE_arDCA_input.csv"
seqdesign_file="data/clustered_sampling/UBE_100_mutants.csv"
EV=pd.read_csv(EV_file)

### process arDCA to get correct mutations 
arDCA=pd.read_csv(arDCA_file)
arDCA=arDCA.set_index("AACombo", drop=False)
arDCA_input=pd.read_csv(arDCA_input)
arDCA_input=arDCA_input.set_index('converted_mutations', drop=False)
arDCA=arDCA.join(arDCA_input)
arDCA=arDCA[['log_fitness', 'arDCA', 'mutations']]
arDCA=arDCA.set_index("mutations", drop=False)

### 'arDCA', 'AACombo'
seqdesign=pd.read_csv(seqdesign_file)
seqdesign=seqdesign.set_index("name", drop=False)
seqdesign=seqdesign.rename(columns={'mean': 'seqdesign'})
### 'name', 'mean'
df=EV[['log_fitness', 'mutations', 'plmc']]
df = df.rename(columns={'log_fitness': 'Fitness', 'mutations': 'AACombo'})
df=df.set_index("AACombo", drop=False)
df=df.join(arDCA[["arDCA",'log_fitness']])
df=df.join(seqdesign[["name", "seqdesign"]])

features=np.load(encoding_file)
ComboToIndex=pickle.load(open(ComboToIndex_file, "rb"))
AACombo=list(ComboToIndex.keys())
if len(features.shape) == 3:
    features = np.reshape(features, [features.shape[0], features.shape[1] * features.shape[2]])

distortions = []
K = range(1,30, 2)
for k in K:
    kmeanModel = KMeans(n_clusters=k)
    kmeanModel.fit(features)
    distortions.append(kmeanModel.inertia_)

plt.figure(figsize=(16,8))
plt.plot(K, distortions, 'bx-')
#plt.axvline(x=10, color='k', linestyle='--')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('The Elbow Method showing the optimal k')
plt.show()
### set cluster to 30
kmeans = KMeans(n_clusters=27).fit(features)
cluster_labels = kmeans.labels_
clusters=pd.DataFrame(list(zip(AACombo, kmeans.labels_)), columns=["AACombo", "cluster_id"])
clusters=clusters.set_index("AACombo", drop=False)
df=df.join(clusters["cluster_id"])
df["Fitness_n"]=(df["Fitness"].values - min(df["Fitness"].values))/ (max(df["Fitness"].values)-min(df["Fitness"].values))
df[['AACombo', 'Fitness', 'Fitness_n', 'plmc', 'arDCA','seqdesign', 'cluster_id']].\
to_csv("data/clustered_sampling/UBE_zeroshot.csv")

### PhoQ
encoding_file="PhoQ_Georgiev_normalized.npy"
ComboToIndex_file="ComboToIndex_PhoQ.pkl"
EV_file="data/clustered_sampling/PhoQ_fitness_EV.csv"
arDCA_file="data/clustered_sampling/PhoQ_arDCA.csv"
seqdesign_file="data/clustered_sampling/PhoQ_mutant_100_1.csv"
EV=pd.read_csv(EV_file)

### process arDCA to get correct mutations 
arDCA=pd.read_csv(arDCA_file)
arDCA=arDCA.set_index("AACombo", drop=False)
### 'arDCA', 'AACombo'
seqdesign=pd.read_csv(seqdesign_file)
seqdesign=seqdesign.set_index("name", drop=False)
seqdesign=seqdesign.rename(columns={'mean': 'seqdesign'})
### 'name', 'mean'
df=EV[['Variants', 'Fitness', 'plmc']]
df = df.rename(columns={'Variants': 'AACombo'})
df=df.set_index("AACombo", drop=False)
df=df.join(arDCA[["arDCA"]])
df=df.join(seqdesign[["name", "seqdesign"]])

features=np.load(encoding_file)
ComboToIndex=pickle.load(open(ComboToIndex_file, "rb"))
AACombo=list(ComboToIndex.keys())
if len(features.shape) == 3:
    features = np.reshape(features, [features.shape[0], features.shape[1] * features.shape[2]])

distortions = []
K = range(1,30, 2)
for k in K:
    kmeanModel = KMeans(n_clusters=k)
    kmeanModel.fit(features)
    distortions.append(kmeanModel.inertia_)

plt.figure(figsize=(16,8))
plt.plot(K, distortions, 'bx-')
#plt.axvline(x=10, color='k', linestyle='--')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('The Elbow Method showing the optimal k')
plt.show()
### set cluster to 30
kmeans = KMeans(n_clusters=19).fit(features)
cluster_labels = kmeans.labels_
clusters=pd.DataFrame(list(zip(AACombo, kmeans.labels_)), columns=["AACombo", "cluster_id"])
clusters=clusters.set_index("AACombo", drop=False)
df=df.join(clusters["cluster_id"])
df["Fitness_n"]=(df["Fitness"].values - min(df["Fitness"].values))/ (max(df["Fitness"].values)-min(df["Fitness"].values))
df[['AACombo', 'Fitness', 'Fitness_n', 'plmc', 'arDCA','seqdesign', 'cluster_id']].\
to_csv("data/clustered_sampling/PhoQ_zeroshot.csv")

### SaCas9_WED
encoding_file="SaCas9_888889909_georgiev_UnNormalized.npy"
ComboToIndex_file="SaCas9_888889909_georgiev_ComboToIndex.pkl"
EV_file="data/clustered_sampling/SaCas9_WED_EV.csv"
arDCA_file="data/clustered_sampling/SaCas9_N888_A889_L909_arDCA.csv"
seqdesign_file="data/clustered_sampling/SaCas9_888889909_mutant_100.csv"
EV=pd.read_csv(EV_file)
### process arDCA to get correct mutations 
arDCA=pd.read_csv(arDCA_file)
arDCA=arDCA.set_index("AACombo", drop=False)
### 'arDCA', 'AACombo'
seqdesign=pd.read_csv(seqdesign_file)
seqdesign=seqdesign.set_index("name", drop=False)
seqdesign=seqdesign.rename(columns={'mean': 'seqdesign'})
### 'name', 'mean'
df=EV
df = df.rename(columns={'Variants': 'AACombo'})
df=df.set_index("AACombo", drop=False)
df=df.join(arDCA[["arDCA"]])
df=df.join(seqdesign[["name", "seqdesign"]])

features=np.load(encoding_file)
ComboToIndex=pickle.load(open(ComboToIndex_file, "rb"))
AACombo=list(ComboToIndex.keys())
if len(features.shape) == 3:
    features = np.reshape(features, [features.shape[0], features.shape[1] * features.shape[2]])

distortions = []
K = range(1,30, 2)
for k in K:
    kmeanModel = KMeans(n_clusters=k)
    kmeanModel.fit(features)
    distortions.append(kmeanModel.inertia_)

plt.figure(figsize=(16,8))
plt.plot(K, distortions, 'bx-')
#plt.axvline(x=10, color='k', linestyle='--')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('The Elbow Method showing the optimal k')
plt.show()
### set cluster to 25
kmeans = KMeans(n_clusters=25).fit(features)
cluster_labels = kmeans.labels_
clusters=pd.DataFrame(list(zip(AACombo, kmeans.labels_)), columns=["AACombo", "cluster_id"])
clusters=clusters.set_index("AACombo", drop=False)
df=df.join(clusters["cluster_id"])
df[['AACombo', 'sgRNA', 'Fitness_n', 'plmc', 'arDCA','seqdesign', 'cluster_id']].\
to_csv("SaCas9_WED_zeroshot.csv")


