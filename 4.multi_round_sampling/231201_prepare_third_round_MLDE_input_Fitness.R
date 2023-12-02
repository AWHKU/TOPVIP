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

outdir<-"data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/3rd_round_MLDE_input_Fitness_files/"
### SpCas9 ###
Fit_norm<-read_csv("data/MLDE_train_test_datasets/SpCas9_fitness_norm.csv")
Fit_norm<- Fit_norm %>% gather(sg, Fitness_n, 2:3)
## input data in Supplementary table S1
Fit<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SpCas9_Fitness.csv")
Fit<-left_join(Fit, Fit_norm)
test_sample<-read_csv("data/MLDE_train_test_datasets/SpCas9_test_sample.data")

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SpCas9_Sg5_arDCA_CLADETOP_georgiev_p2/20230618-222252/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="Sg5") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}

o<- bind_rows(o, Fit%>% filter(sg=="Sg5") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SpCas9_Sg5_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SpCas9_Sg5_EVarDCA_CLADE2TOP_georgiev_p2/20230627-141131/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="Sg5") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}

o<- bind_rows(o, Fit%>% filter(sg=="Sg5") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SpCas9_Sg5_EVarDCA_CLADE2TOPCLADE.csv"))


t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SpCas9_Sg5_EVmutation_CLADETOP_georgiev_p2/20230618-222411/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="Sg5") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}

o<- bind_rows(o, Fit%>% filter(sg=="Sg5") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SpCas9_Sg5_EVmutation_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SpCas9_Sg8_arDCA_CLADETOP_georgiev_p2/20230618-222516/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="Sg8") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}

o<- bind_rows(o, Fit%>% filter(sg=="Sg8") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SpCas9_Sg8_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SpCas9_Sg8_EVarDCA_CLADE2TOP_georgiev_p2/20230627-141238/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="Sg8") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}

o<- bind_rows(o, Fit%>% filter(sg=="Sg8") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SpCas9_Sg8_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SpCas9_Sg8_EVmutation_CLADETOP_georgiev_p2/20230618-222637/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="Sg8") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}

o<- bind_rows(o, Fit%>% filter(sg=="Sg8") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SpCas9_Sg8_EVmutation_CLADETOPCLADE.csv"))

### SaCas9 1296 ###
Fit_norm<-read_csv("data/MLDE_train_test_datasets/SaCas9_887888887985986988989991_Fitness_norm.csv")
Fit_norm<- Fit_norm %>% gather(sg, Fitness_n, 2:4)
Fit<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/Sa1296_Fitness.csv")
Fit<-left_join(Fit, Fit_norm)
test_sample<-read_csv("data/MLDE_train_test_datasets/SaCas9_887888889985986988989991_fixed_test_set.txt")

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas91296_ACAAGT_arDCA_CLADETOP_georgiev_p2/20230618-223409/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="ACAAGT") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="ACAAGT") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SaCas91296_ACAAGT_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas91296_ACAAGT_EVarDCA_CLADE2TOP_georgiev_p2/20230627-142101/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="ACAAGT") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="ACAAGT") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SaCas91296_ACAAGT_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas91296_GGTGGT_arDCA_CLADETOP_georgiev_p2/20230618-223647/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="GGTGGT") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="GGTGGT") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SaCas91296_GGTGGT_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas91296_GGTGGT_EVarDCA_CLADE2TOP_georgiev_p2/20230627-142224/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="GGTGGT") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="GGTGGT") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SaCas91296_GGTGGT_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas91296_TGGAGT_arDCA_CLADETOP_georgiev_p2/20230618-223925/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="TGGAGT") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="TGGAGT") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SaCas91296_TGGAGT_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas91296_TGGAGT_EVarDCA_CLADE2TOP_georgiev_p2/20230627-142329/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="TGGAGT") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="TGGAGT") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SaCas91296_TGGAGT_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas91296_ACAAGT_EVmutation_CLADETOP_georgiev_p2/20230618-223523/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="ACAAGT") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="ACAAGT") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SaCas91296_ACAAGT_EVmutation_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas91296_GGTGGT_EVmutation_CLADETOP_georgiev_p2/20230618-223817/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="GGTGGT") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="GGTGGT") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SaCas91296_GGTGGT_EVmutation_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas91296_TGGAGT_EVmutation_CLADETOP_georgiev_p2/20230618-224028/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="TGGAGT") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="TGGAGT") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o, paste0(outdir, "SaCas91296_TGGAGT_EVmutation_CLADETOPCLADE.csv"))

### SaCas9 8000 ###
Fit_norm<-read_csv("data/MLDE_train_test_datasets/SaCas9_8888889909_Fitness_norm.csv")
Fit_norm<- Fit_norm %>% gather(sg, Fitness_n, 2:5)
Fit<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/Sa8000_Fitness.csv")
Fit<-left_join(Fit, Fit_norm)
test_sample<-read_csv("data/MLDE_train_test_datasets/Sa8000_test_sample2.data")

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp56libD18_arDCA_CLADETOP_georgiev_p2/20230627-141455/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp56libD18") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp56libD18") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp56libD18_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp58libD21_arDCA_CLADETOP_georgiev_p2/20230618-223224/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp58libD21") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp58libD21") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp58libD21_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp55libD18_EVmutation_CLADETOP_georgiev_p2/20230618-224327/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp55libD18") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp55libD18") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp55libD18_EVmutation_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp56libD18_EVmutation_CLADETOP_georgiev_p2/20230627-141653/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp56libD18") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp56libD18") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp56libD18_EVmutation_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp58libD14_EVmutation_CLADETOP_georgiev_p2/20230618-224540/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp58libD14") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp58libD14") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp58libD14_EVmutation_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp58libD21_EVmutation_CLADETOP_georgiev_p2/20230618-224653/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp58libD21") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp58libD21") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp58libD21_EVmutation_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp55libD18_EVarDCA_CLADE2TOP_georgiev_p2/20230627-141345/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp55libD18") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp55libD18") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp55libD18_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp56libD18_EVarDCA_CLADE2TOP_georgiev_p2/20230627-141553/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp56libD18") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp56libD18") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp56libD18_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp58libD14_EVarDCA_CLADE2TOP_georgiev_p2/20230627-141748/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp58libD14") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp58libD14") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp58libD14_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/SaCas98000_DTp58libD21_EVarDCA_CLADE2TOP_georgiev_p2/20230627-141851/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% filter(sg=="DTp58libD21") %>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(sg=="DTp58libD21") %>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "SaCas98000_DTp58libD21_EVarDCA_CLADE2TOPCLADE.csv"))

### PhoQ ###
Fit<-read_csv("data/MLDE_train_test_datasets/PhoQ_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
test_sample<-read_csv("data/MLDE_train_test_datasets/PhoQ_test_sample.data")
colnames(test_sample)<-c("AACombo", "Fitness", "Fitness_n")

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/PhoQ_arDCA_CLADETOP_georgiev_p2/20230618-221638/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "PhoQ_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/PhoQ_EVarDCA_CLADE2TOP_georgiev_p2/20230627-140849/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "PhoQ_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/PhoQ_EVmutation_CLADETOP_georgiev_p2/20230618-221809/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "PhoQ_EVmutation_CLADETOPCLADE.csv"))

### GFP ###
Fit<-read_csv("data/MLDE_train_test_datasets/GFP_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
test_sample<-read_csv("data/MLDE_train_test_datasetsGFP_test_sample.data")

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/GFP_arDCA_CLADETOP_georgiev_p2/20230618-215843/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "GFP_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/GFP_EVarDCA_CLADE2TOP_georgiev_p2/20230627-140055/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "GFP_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/GFP_EVmutation_CLADETOP_georgiev_p2/20230618-220039/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "GFP_EVmutation_CLADETOPCLADE.csv"))

### PABP ###
Fit<-read_csv("data/MLDE_train_test_datasets/PABP_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
test_sample<-read_csv("data/MLDE_train_test_datasets/PABP_test_sample.data")

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/PABP_arDCA_CLADETOP_georgiev_p2/20230618-220740/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "PABP_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/PABP_EVarDCA_CLADE2TOP_georgiev_p2/20230627-140327/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "PABP_EVarDCA_CLADE2TOPCLADE.csv"))
t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/PABP_EVmutation_CLADETOP_georgiev_p2/20230618-220926/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "PABP_EVmutation_CLADETOPCLADE.csv"))

### UBE ###
Fit<-read_csv("data/MLDE_train_test_datasets/UBE_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
test_sample<-read_csv("data/MLDE_train_test_datasets/UBE_test_sample.data")

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/UBE_arDCA_CLADETOP_georgiev_p2/20230618-221159/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "UBE_arDCA_CLADETOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/UBE_EVarDCA_CLADE2TOP_georgiev_p2/20230627-140545/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "UBE_EVarDCA_CLADE2TOPCLADE.csv"))

t<-read_csv("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark/UBE_EVmutation_CLADETOP_georgiev_p2/20230618-221404/PredictedFitness.csv")
training_variants<- t %>% filter(`InTrainingData?`=="YES") %>% select(AACombo) %>% unlist()
t<-t %>% filter(`InTrainingData?`=="NO") %>% filter(!AACombo %in% test_sample)
t<-left_join(t, Fit%>% select(AACombo, cluster_id, Fitness_n))
t<-t %>% filter(!is.na(Fitness_n))
Prob<-t %>% group_by(cluster_id) %>% summarise(m=mean(PredictedFitness, na.rm=T)) %>% mutate(p=m/sum(m)) %>% select(cluster_id, p)
Prob$p[Prob$p<0]<-min(Prob$p[Prob$p>0])
Prob<-as_tibble(sample(Prob$cluster_id, 12, replace = T, prob=Prob$p) %>% table())
colnames(Prob)<-c("cluster_id", "n")
o<-tibble()
for (i in 1:nrow(Prob)){
  cluster=as.numeric(Prob[i, 1])
  n=unlist(Prob[i, 2])
  o<-bind_rows(o,
               t %>% filter(cluster_id==cluster) %>% arrange(-PredictedFitness) %>% head(n=n) %>% select(AACombo, Fitness_n))
}
o<- bind_rows(o, Fit%>% filter(AACombo %in% training_variants) %>% select(AACombo, cluster_id, Fitness_n))
o<- o %>% ungroup() %>% select(AACombo, Fitness_n)
colnames(o)<-c("AACombo", "Fitness")
write_csv(o %>% arrange(-Fitness) %>% head(n=36) , paste0(outdir, "UBE_EVmutation_CLADETOPCLADE.csv"))
