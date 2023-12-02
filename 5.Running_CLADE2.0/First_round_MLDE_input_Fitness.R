library(seqinr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggrepel)
library(ggsci)
theme_publication<-theme_classic()+
  theme(text = element_text(size=10, family = "Arial"),
        panel.background = element_rect(colour="black", linewidth=0.5),
        legend.position = "bottom",
        strip.background = element_blank())

outdir<-"data/execute_MLDE_benchmark/1st_round_MLDE_input_Fitness_files"
### SpCas9 ###
### CLADE2.0 in put ###
Fit_norm<-read_csv("data/MLDE_train_test_datasets/SpCas9_fitness_norm.csv")
Fit_norm<- Fit_norm %>% gather(sg, Fitness_n, 2:3)
## input data in Supplementary table S1
Fit<-read_csv("data/execute_MLDE_benchmark/SpCas9_Fitness.csv")

Fit<-left_join(Fit, Fit_norm)
test_sample<-read_csv("data/MLDE_train_test_datasets/SpCas9_test_sample.data")
### EV
t<-read_csv("data/execute_MLDE_benchmark/SpCas9_EV_InputValidationData.csv")

### 3 x 48 random for augmented representation 
o<-Fit %>% filter(sg=="Sg5") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>%  select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "SpCas9_Sg5_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SpCas9_Sg5_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SpCas9_Sg5_random3_Fitness.csv"))
o<-Fit %>% filter(sg=="Sg8") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "SpCas9_Sg8_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SpCas9_Sg8_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SpCas9_Sg8_random3_Fitness.csv"))
# 1 x 48 top seqdesign
o<-Fit %>% filter(sg=="Sg5") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SpCas9_Sg5_seqdesign_Fitness.csv"))
o<-t %>% filter(sg=="Sg8") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SpCas9_Sg8_seqdesign_Fitness.csv"))
### 12 x EVmutation clade ## use cluster ID in Fit
WT<-"RQKETQKA"
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(EvMutation, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="Sg5") %>% filter(cluster_id==cluster) %>% arrange(desc(EvMutation)) %>% head(n=n) %>% select(AACombo, Fitness_n))
  }
if (!WT %in% o$AACombo){
 o<-bind_rows(o, 
              t %>% filter(AACombo ==WT, sg=="Sg5") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SpCas9_Sg5_CLADEEVmutation_Fitness.csv"))

o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="Sg8") %>% filter(cluster_id==cluster) %>% arrange(desc(EvMutation)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="Sg8") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SpCas9_Sg8_CLADEEVmutation_Fitness.csv"))

### 12 x arDCA clade 
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(-arDCA, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
  t %>% filter(sg=="Sg5") %>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="Sg5") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SpCas9_Sg5_CLADEarDCA_Fitness.csv"))

o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="Sg8") %>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o$AACombo){
  o<-bind_rows(o[1:12, ], 
               t %>% filter(AACombo ==WT, sg=="Sg8") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}

colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SpCas9_Sg8_CLADEarDCA_Fitness.csv"))
### SaCas9_1296 ###
### CLADE2.0 in put ###
Fit_norm<-read_csv("data/MLDE_train_test_datasets/SaCas9_887888887985986988989991_Fitness_norm.csv")
Fit_norm<- Fit_norm %>% gather(sg, Fitness_n, 2:4)
Fit<-read_csv("data/execute_MLDE_benchmark/Sa1296_Fitness.csv")
Fit<-left_join(Fit, Fit_norm)
test_sample<-read_csv("data/MLDE_train_test_datasets/SaCas9_887888889985986988989991_fixed_test_set.txt")
### EV
t<-read_csv("data/execute_MLDE_benchmark/SaCas9_1296_EV_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness.y))

o<-t %>% filter(sg=="ACAAGT") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_ACAAGT_CLADE2EV_Fitness.csv"))
o<-t %>% filter(sg=="GGTGGT") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_GGTGGT_CLADE2EV_Fitness.csv"))
o<-t %>% filter(sg=="TGGAGT") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_TGGAGT_CLADE2EV_Fitness.csv"))
### EV + seqdesign 
t<-read_csv("data/execute_MLDE_benchmark/SaCas9_1296_EV+seqdesign_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness.y))

o<-t %>% filter(sg=="ACAAGT") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_ACAAGT_CLADE2EV+seqdesign__Fitness.csv"))
o<-t %>% filter(sg=="GGTGGT") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_GGTGGT_CLADE2EV+seqdesign__Fitness.csv"))
o<-t %>% filter(sg=="TGGAGT") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_TGGAGT_CLADE2EV+seqdesign__Fitness.csv"))

### 3 x 48 random for augmented representation 
o<-Fit %>% filter(sg=="ACAAGT") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_ACAAGT_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_ACAAGT_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_ACAAGT_random3_Fitness.csv"))
o<-Fit %>% filter(sg=="GGTGGT") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_GGTGGT_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_GGTGGT_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_GGTGGT_random3_Fitness.csv"))
o<-Fit %>% filter(sg=="TGGAGT") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_TGGAGT_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_TGGAGT_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_TGGAGT_random3_Fitness.csv"))

# 1 x 48 top seqdesign
o<-Fit %>% filter(sg=="ACAAGT") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_ACAAGT_seqdesign_Fitness.csv"))
o<-Fit %>% filter(sg=="GGTGGT") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_GGTGGT_seqdesign_Fitness.csv"))
o<-Fit %>% filter(sg=="TGGAGT") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_TGGAGT_seqdesign_Fitness.csv"))
### 12 x EVmutation clade ## use cluster ID in Fit
WT<-"LNANNLLR"
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(EvMutation, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="ACAAGT") %>% filter(cluster_id==cluster) %>% arrange(desc(EvMutation)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="ACAAGT") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_ACAAGT_CLADEEVmutation_Fitness.csv"))

o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="GGTGGT") %>% filter(cluster_id==cluster) %>% arrange(desc(EvMutation)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="GGTGGT") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_GGTGGT_CLADEEVmutation_Fitness.csv"))

o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="TGGAGT") %>% filter(cluster_id==cluster) %>% arrange(desc(EvMutation)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="TGGAGT") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_TGGAGT_CLADEEVmutation_Fitness.csv"))


### 12 x arDCA clade 
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(-arDCA, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="ACAAGT") %>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o$AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="ACAAGT") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_ACAAGT_CLADEEarDCA_Fitness.csv"))

o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="GGTGGT") %>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o$AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="GGTGGT") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_GGTGGT_CLADEEarDCA_Fitness.csv"))

o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="TGGAGT") %>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o$AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="TGGAGT") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_TGGAGT_CLADEarDCA_Fitness.csv"))

### SaCas9_888889909
### CLADE2.0 in put ###
Fit_norm<-read_csv("data/MLDE_train_test_datasets/SaCas9_8888889909_Fitness_norm.csv")
Fit_norm<- Fit_norm %>% gather(sg, Fitness_n, 2:5)
Fit<-read_csv("data/execute_MLDE_benchmark/Sa8000_Fitness.csv")
Fit<-left_join(Fit, Fit_norm)
test_sample<-read_csv("data/MLDE_train_test_datasets/Sa8000_test_sample2.data")
### EV
t<-read_csv("data/execute_MLDE_benchmark/SaCas9_8000_EV_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness.y))

o<-t %>% filter(sg=="DTp56libD18") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_DTp56libD18_CLADE2EV_Fitness.csv"))
o<-t %>% filter(sg=="DTp58libD21") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_DTp58libD21_CLADE2EV_Fitness.csv"))

### EV + seqdesign 
t<-read_csv("data/execute_MLDE_benchmark/SaCas9_8000_EV+seqdesign_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness.y))

o<-t %>% filter(sg=="DTp56libD18") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_DTp56libD18_CLADE2EV+seqdesign_Fitness.csv"))
o<-t %>% filter(sg=="DTp58libD21") %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_DTp58libD21_CLADE2EV+seqdesign_Fitness.csv"))

### 3 x 48 random for augmented representation 
o<-Fit %>% filter(sg=="DTp56libD18") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_DTp56libD18_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_DTp56libD18_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_DTp56libD18_random3_Fitness.csv"))
o<-Fit %>% filter(sg=="DTp58libD21") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_DTp58libD21_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_DTp58libD21_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "SaCas9_DTp58libD21_random3_Fitness.csv"))

# 1 x 48 top seqdesign
o<-Fit %>% filter(sg=="DTp56libD18") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_DTp56libD18_seqdesign_Fitness.csv"))
o<-Fit %>% filter(sg=="DTp58libD21") %>% filter(!is.na(Fitness_n)) %>% filter(! AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "SaCas9_DTp58libD21_seqdesign_Fitness.csv"))

### 12 x EVmutation clade ## use cluster ID in Fit
WT<-"NAL"
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(EvMutation, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="DTp56libD18") %>% filter(cluster_id==cluster) %>% arrange(desc(EvMutation)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="DTp56libD18") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_DTp56libD18_CLADEEVmutation_Fitness.csv"))
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="DTp58libD21") %>% filter(cluster_id==cluster) %>% arrange(desc(EvMutation)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="DTp58libD21") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_DTp58libD21_CLADEEVmutation_Fitness.csv"))

### 12 x arDCA clade 
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(-arDCA, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="DTp56libD18") %>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="DTp56libD18") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_DTp56libD18_CLADEEarDCA_Fitness.csv"))
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(sg=="DTp58libD21") %>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT, sg=="DTp58libD21") %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "SaCas9_DTp58libD21_CLADEarDCA_Fitness.csv"))


### PhoQ ### 
### CLADE2.0 in put ###
Fit<-read_csv("data/MLDE_train_test_datasets/PhoQ_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
test_sample<-read_csv("data/MLDE_train_test_datasets/PhoQ_test_sample.data")
colnames(test_sample)<-c("AACombo", "Fitness", "Fitness_n")
### EV
t<-read_csv("data/execute_MLDE_benchmark/PhoQ_EV_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness_n))

o<-t %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "PhoQ_CLADE2EV_Fitness.csv"))
### EV + seqdesign 
t<-read_csv("data/execute_MLDE_benchmark/PhoQ_EV+seqdesign_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness_n))

o<-t %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "PhoQ_CLADE2EV+seqdesign_Fitness.csv"))

### 3 x 48 random for augmented representation 
o<-Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "PhoQ_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "PhoQ_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "PhoQ_random3_Fitness.csv"))

# 1 x 48 top seqdesign
o<-Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "PhoQ_seqdesign_Fitness.csv"))

### 12 x EVmutation clade ## use cluster ID in Fit
WT<-"AVST"
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(plmc, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0], na.rm = T)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t  %>% filter(cluster_id==cluster) %>% arrange(desc(plmc)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t  %>% filter(AACombo ==WT) %>%  select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "PhoQ_CLADEEVmutation_Fitness.csv"))

### 12 x arDCA clade 
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(-arDCA, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0], na.rm = T)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t%>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT) %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "PhoQ_CLADEEarDCA_Fitness.csv"))

### GFP ###
### CLADE2.0 in put ###
Fit<-read_csv("data/MLDE_train_test_datasets/GFP_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
test_sample<-read_csv("data/MLDE_train_test_datasetsGFP_test_sample.data")
### EV
t<-read_csv("data/execute_MLDE_benchmark/GFP_EV_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness_n))

o<-t %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "GFP_CLADE2EV_Fitness.csv"))
### EV + seqdesign 
t<-read_csv("data/execute_MLDE_benchmark/GFP_EV+seqdesign_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness_n))

o<-t %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "GFP_CLADE2EV+seqdesign_Fitness.csv"))

### 3 x 48 random for augmented representation 
o<-Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "GFP_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "GFP_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "GFP_random3_Fitness.csv"))

# 1 x 48 top seqdesign
o<-Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "GFP_seqdesign_Fitness.csv"))

### 12 x EVmutation clade ## use cluster ID in Fit
WT<-"WT"
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(plmc, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0], na.rm = T)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t  %>% filter(cluster_id==cluster) %>% arrange(desc(plmc)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t  %>% filter(AACombo ==WT) %>%  select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "GFP_CLADEEVmutation_Fitness.csv"))

### 12 x arDCA clade 
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(-arDCA, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0], na.rm = T)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t%>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o $AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT) %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "GFP_CLADEEarDCA_Fitness.csv"))

### PABP ###
### CLADE2.0 in put ###
Fit<-read_csv("data/MLDE_train_test_datasets/PABP_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
test_sample<-read_csv("data/MLDE_train_test_datasets/PABP_test_sample.data")
### EV
t<-read_csv("data/execute_MLDE_benchmark/PABP_EV_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness_n))

o<-t %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "PABP_CLADE2EV_Fitness.csv"))
### EV + seqdesign 
t<-read_csv("data/execute_MLDE_benchmark/PABP_EV+seqdesign_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness_n))

o<-t %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "PABP_CLADE2EV+seqdesign_Fitness.csv"))

### 3 x 48 random for augmented representation 
o<-Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "PABP_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "PABP_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "PABP_random3_Fitness.csv"))

# 1 x 48 top seqdesign
o<-Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "PABP_seqdesign_Fitness.csv"))

### 12 x EVmutation clade ## use cluster ID in Fit
WT<-"WT"
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(plmc, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0], na.rm = T)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t  %>% filter(cluster_id==cluster) %>% arrange(desc(plmc)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o$AACombo){
  o<-bind_rows(o, 
               t  %>% filter(AACombo ==WT) %>%  select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "PABP_CLADEEVmutation_Fitness.csv"))

### 12 x arDCA clade 
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(-arDCA, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0], na.rm = T)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t%>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o$AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT) %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "PABP_CLADEEarDCA_Fitness.csv"))

### UBE ###
### CLADE2.0 in put ###
Fit<-read_csv("data/MLDE_train_test_datasets/UBE_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
test_sample<-read_csv("data/MLDE_train_test_datasets/UBE_test_sample.data")
### EV
t<-read_csv("data/execute_MLDE_benchmark/UBE_EV_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness_n))

o<-t %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "UBE_CLADE2EV_Fitness.csv"))
### EV + seqdesign 
t<-read_csv("data/execute_MLDE_benchmark/UBE_EV_InputValidationData.csv")
### remove test
t<-t %>% mutate(index=seq(1, nrow(t), 1))
t<- t %>% filter(! AACombo %in% test_sample$AACombo)
t<- left_join(t, Fit, by="AACombo")
t<- t %>% filter(!is.na(Fitness_n))

o<-t %>% arrange(index) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "UBE_CLADE2EV+seqdesign_Fitness.csv"))

### 3 x 48 random for augmented representation 
o<-Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo) %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo",   "Fitness")
write_csv(sample_n(o, 48), paste0(outdir, "UBE_random1_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "UBE_random2_Fitness.csv"))
write_csv(sample_n(o, 48), paste0(outdir, "UBE_random3_Fitness.csv"))

# 1 x 48 top seqdesign
o<-Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo) %>% arrange(desc(seqdesign)) %>% select(AACombo, Fitness_n)
o<-o[1:48, ]
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o, paste0(outdir, "UBE_seqdesign_Fitness.csv"))

### 12 x EVmutation clade ## use cluster ID in Fit
WT<-"WT"
### UBE has no WT fitness
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(plmc, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0], na.rm = T)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t  %>% filter(cluster_id==cluster) %>% arrange(desc(plmc)) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o$AACombo){
  o<-bind_rows(o, 
               t  %>% filter(AACombo ==WT) %>%  select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "UBE_CLADEEVmutation_Fitness.csv"))

### 12 x arDCA clade 
Prob<-Fit %>% group_by(cluster_id) %>% summarise(m=mean(-arDCA, na.rm=T)) 
Prob$m[is.na(Prob$m)]<-min(Prob$m, na.rm = T)
Prob<-Prob%>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0], na.rm = T)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
t<- Fit %>% filter(!is.na(Fitness_n)) %>% filter(!AACombo %in% test_sample$AACombo)
### select top variants per  cluster ###
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t%>% filter(cluster_id==cluster) %>% arrange(arDCA) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
if (!WT %in% o$AACombo){
  o<-bind_rows(o, 
               t %>% filter(AACombo ==WT) %>% select(AACombo, Fitness_n)) %>% arrange(desc(Fitness_n))
}
colnames(o)<-c("AACombo",   "Fitness")
write_csv(o[1:12, ], paste0(outdir, "UBE_CLADEEarDCA_Fitness.csv"))

###  final check, make sure all input fitness files do not have test_sample AACombo ###
test_sample<-read_csv("data/MLDE_train_test_datasets/SpCas9_R661Q695K848E923T924Q926K1003R1060_fixed_test_set.txt")
fit_files<-list.files(path=outdir, pattern="SpCas9")
for (f in fit_files){
  f_<-read_csv(paste0(outdir, f), )
  r<-f_ %>% filter(AACombo %in% test_sample$AACombo) %>% nrow()
  if (r>0){
  print(c( f, " contain test data " , r))
    }
}

test_sample<-read_csv("data/MLDE_train_test_datasets/SaCas9_887888889985986988989991_fixed_test_set.txt")
fit_files<-list.files(path=outdir, pattern="SaCas9")
fit_files<-fit_files[!grepl("SaCas9_DT", fit_files)]
for (f in fit_files){
  f_<-read_csv(paste0(outdir, f), )
  r<-f_ %>% filter(AACombo %in% test_sample$AACombo) %>% nrow()
  if (r>0){
    print(c( f, " contain test data " , r))
  }
}

test_sample<-read_csv("data/MLDE_train_test_datasets/Sa8000_test_sample2.data")
fit_files<-list.files(path=outdir, pattern="SaCas9")
fit_files<-fit_files[grepl("SaCas9_DT", fit_files)]
for (f in fit_files){
  f_<-read_csv(paste0(outdir, f), )
  r<-f_ %>% filter(AACombo %in% test_sample$AACombo) %>% nrow()
  if (r>0){
    print(c( f, " contain test data " , r))
  }
}

test_sample<-read_csv("data/MLDE_train_test_datasets/PhoQ_test_sample.data")
fit_files<-list.files(path=outdir, pattern="PhoQ")
for (f in fit_files){
  f_<-read_csv(paste0(outdir, f), )
  r<-f_ %>% filter(AACombo %in% test_sample$Variants) %>% nrow()
  if (r>0){
    print(c( f, " contain test data " , r))
  }
}

test_sample<-read_csv("data/MLDE_train_test_datasets/evo3000_test_sample.data")
fit_files<-list.files(path=outdir, pattern="evo")
for (f in fit_files){
  f_<-read_csv(paste0(outdir, f), )
  r<-f_ %>% filter(AACombo %in% test_sample$AACombo) %>% nrow()
  if (r>0){
    print(c( f, " contain test data " , r))
  }
}

test_sample<-read_csv("data/MLDE_train_test_datasets/PABP_test_sample.data")
fit_files<-list.files(path=outdir, pattern="PABP")
for (f in fit_files){
  f_<-read_csv(paste0(outdir, f), )
  r<-f_ %>% filter(AACombo %in% test_sample$AACombo) %>% nrow()
  if (r>0){
    print(c( f, " contain test data " , r))
  }
}

test_sample<-read_csv("data/MLDE_train_test_datasets/GFP_test_sample.data")
fit_files<-list.files(path=outdir, pattern="GFP")
for (f in fit_files){
  f_<-read_csv(paste0(outdir, f), )
  r<-f_ %>% filter(AACombo %in% test_sample$AACombo) %>% nrow()
  if (r>0){
    print(c( f, " contain test data " , r))
  }
}


test_sample<-read_csv("data/MLDE_train_test_datasets/UBE_test_sample.data")
fit_files<-list.files(path=outdir, pattern="UBE")
for (f in fit_files){
  f_<-read_csv(paste0(outdir, f), )
  r<-f_ %>% filter(AACombo %in% test_sample$AACombo) %>% nrow()
  if (r>0){
    print(c( f, " contain test data " , r))
  }
}
