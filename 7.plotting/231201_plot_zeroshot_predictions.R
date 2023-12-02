library(ggplot2)
library(tidyverse)
library(ggpmisc)
library(ggseqlogo)
library(readxl)
library(Cairo)
library(stringdist)
library(ggsci)
library(ggrepel)
theme_publication<-theme_classic()+
  theme(text = element_text(size=12, family = "Arial"),
        panel.background = element_rect(colour="black", size=0.5),
        legend.position = "bottom",
        strip.background = element_blank())

### SpCas9
### input file: Supplementary table 1 ###
t<-read_excel(skip=1, "tableS1.xlsx", sheet = 1)
t<-t %>% gather(sg, Fitness, 2:3)
cor_eff<-c()
metrics<-c()
sg<-c()
g<-t%>% filter(sg=="Sg5") 
for (i in 2:6){
  cor_eff<- c(cor_eff, cor.test(g[, "Fitness"] %>% unlist(), g[, i] %>% unlist(), method="spearman")$estimate)
  metrics<-c(metrics, colnames(g)[i])
  sg<-c(sg, "Sg5")
}
g<-g %>%gather(metrics, values, 2:6)
g$values[g$metrics %in% c("EvoEF", "arDCA")]<- -1 * g$values[g$metrics %in% c("EvoEF", "arDCA")]
tab<- left_join(
  left_join(g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
  filter(rank_n>=0.95, rank_p>=0.95) %>% group_by(sg, metrics) %>% summarise(n=n()),
  g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
    filter(rank_n>=0.95) %>% group_by(sg, metrics) %>% summarise(t=n())) %>% mutate(overlap_5=n/t) %>%
  select(sg, metrics, overlap_5),
data.frame(sg, metrics, cor_eff)) %>% mutate(datatype="SpCas9")

cor_eff<-c()
metrics<-c()
sg<-c()
g<-t%>% filter(sg=="Sg8") 
for (i in 2:6){
  cor_eff<- c(cor_eff, cor.test(g[, "Fitness"] %>% unlist(), g[, i] %>% unlist(), method="spearman")$estimate)
  metrics<-c(metrics, colnames(g)[i])
  sg<-c(sg, "Sg8")
}
g<-g %>%gather(metrics, values, 2:6)
g$values[g$metrics %in% c("EvoEF", "arDCA")]<- -1 * g$values[g$metrics %in% c("EvoEF", "arDCA")]
tab<- bind_rows(tab,left_join(
 full_join(g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
              filter(rank_n>=0.95, rank_p>=0.95) %>% group_by(sg, metrics) %>% summarise(n=n()),
            g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
              filter(rank_n>=0.95) %>% group_by(sg, metrics) %>% summarise(t=n())) %>% mutate(overlap_5=n/t) %>%
    select(sg, metrics, overlap_5),
  data.frame(sg, metrics, cor_eff))  %>% mutate(datatype="SpCas9")
)

### look at consensus? 
t<-t %>% gather(predictor, values, 2:6)
t$values<-as.numeric(t$values)
t$values[t$predictor %in% c("EvoEF", "arDCA")]<- -1 * t$values[t$predictor %in% c("EvoEF", "arDCA")]
### try to show the selections in graph
### plot 1) boxplot of the predicted top 50 fitness and plot 2) no. of diff from WT
g<-t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  group_by(AACombo, sg) %>% summarise(c=n())

OptiHF="RAAMVAKR"
WT="RQKETQKR"
eSp="RQAETQAA"
Opti="AQKETQHR"

t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  filter(AACombo %in% c(WT, Opti, OptiHF, eSp)) %>%
  group_by(predictor) %>% summarise(n())

g<-g%>% group_by(sg) %>% mutate(rank_c=rank(-c, ties.method = "min")) 
g<-g %>% mutate(WT_d=stringdist(AACombo, WT))
g%>%filter(rank_c<=10) 
Opti %in% g$AACombo[g$rank_c<=10]  
g<- t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  mutate(WT_d=stringdist(AACombo, WT))
g$predictor[grepl("esm", g$predictor)]<-"esm_msa1_t12_100M_UR50S"
dist_plot_1<-g %>% ggplot(aes(x=predictor, y=WT_d))+ geom_boxplot()+
  theme_publication+
  labs(y="AA distance from WT", title="SpCas9")+
  coord_flip()

### SpG
t<-read_excel(skip=1, "tableS1.xlsx", sheet = 4)
t<-  t %>% gather(sg, Fitness, 2:5)
t<-t [, c( "AACombo", "EvoEF","EVmutation", "Prot_bert-Naive", "MSA-Naive", "arDCA","sg","Fitness")]
for (s in unique(t$sg)){
  cor_eff<-c()
  metrics<-c()
  sg<-c()
  g<-t%>% filter(sg==s) 
  for (i in 2:6){
    cor_eff<- c(cor_eff, cor.test(g[, "Fitness"] %>% unlist(), g[, i] %>% unlist(), method="spearman")$estimate)
    metrics<-c(metrics, colnames(g)[i])
    sg<-c(sg, s)
  }
  g<-g %>%gather(metrics, values, 2:6)
  g$values[g$metrics %in% c("EvoEF", "arDCA")]<- -1 * g$values[g$metrics %in% c("EvoEF", "arDCA")]
  b<- right_join(
    left_join(g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
                filter(rank_n>=0.95, rank_p>=0.95) %>% group_by(sg, metrics) %>% summarise(n=n()),
              g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
                filter(rank_n>=0.95) %>% group_by(sg, metrics) %>% summarise(t=n())) %>% mutate(overlap_5=n/t) %>%
      select(sg, metrics, overlap_5),
    data.frame(sg, metrics, cor_eff))  %>% mutate(datatype="SpG")
  tab<-bind_rows(tab, b)
}
tab<-tab[!duplicated(tab), ]
  
### look at consensus? 
t<-t %>% gather(predictor, values, 2:6)
t$values<-as.numeric(t$values)
t$values[t$predictor %in% c("EvoEF", "arDCA")]<- -1 * t$values[t$predictor %in% c("EvoEF", "arDCA")]
### try to show the selections in graph
### plot 1) boxplot of the predicted top 50 fitness and plot 2) no. of diff from WT
g<-t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  group_by(AACombo, sg) %>% summarise(c=n())

g<-g%>% group_by(sg) %>% mutate(rank_c=rank(-c, ties.method = "min")) 
g<-g %>% mutate(WT_d=stringdist(AACombo, WT))
g%>%filter(rank_c<=10) 

## SaCas9_WED-PI
t<-read_excel(skip=1, "tableS1.xlsx", sheet = 2)
t<-t %>% gather(sg, Fitness, 2:4)
for (s in unique(t$sg)){
cor_eff<-c()
metrics<-c()
sg<-c()
g<-t%>% filter(sg==s) 
for (i in 2:6){
  cor_eff<- c(cor_eff, cor.test(g[, "Fitness"] %>% unlist(), g[, i] %>% unlist(), method="spearman")$estimate)
  metrics<-c(metrics, colnames(g)[i])
  sg<-c(sg, s)
}

g<-g %>%gather(metrics, values, 2:6)
g$values[g$metrics %in% c("EvoEF", "arDCA")]<- -1 * g$values[g$metrics %in% c("EvoEF", "arDCA")]
b<- right_join(
  left_join(g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
              filter(rank_n>=0.95, rank_p>=0.95) %>% group_by(sg, metrics) %>% summarise(n=n()),
            g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
              filter(rank_n>=0.95) %>% group_by(sg, metrics) %>% summarise(t=n())) %>% mutate(overlap_5=n/t) %>%
    select(sg, metrics, overlap_5),
  data.frame(sg, metrics, cor_eff)) %>% mutate(datatype="SaCas9_WED-PI")
tab<-bind_rows(tab, b)
}

### look at consensus? 
t<-t %>% gather(predictor, values, 2:6)
t$values<-as.numeric(t$values)
t$values[t$predictor %in% c("EvoEF", "arDCA")]<- -1 * t$values[t$predictor %in% c("EvoEF", "arDCA")]
### try to show the selections in graph
### plot 1) boxplot of the predicted top 50 fitness and plot 2) no. of diff from WT
g<-t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  group_by(AACombo, sg) %>% summarise(c=n())

N888Q="LQANNLLR"
A889S="LNSNNLLR"
N888Q_A889S="LQSNNLLR"
WT="LNANNLLR"

g<-g%>% group_by(sg) %>% mutate(rank_c=rank(-c, ties.method = "min")) 
g<-g %>% mutate(WT_d=stringdist(AACombo, WT))

g<- t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  mutate(WT_d=stringdist(AACombo, WT))

g$predictor[grepl("esm", g$predictor)]<-"esm_msa1_t12_100M_UR50S"
dist_plot_2<-g %>% ggplot(aes(x=predictor, y=WT_d))+ geom_boxplot()+
  theme_publication+
  labs(y="AA distance from WT", title="SaCas9_WED-PI")+
  coord_flip()

WT %in% g$AACombo
N888Q %in% g$AACombo
dist_plot_1+dist_plot_2+plot_layout(ncol = 2)


### SAV SaCas9
t<-read_excel(skip=1, "tableS1.xlsx", sheet = 3)
t<-t %>% gather(sg, Fitness, 2:3)

for (s in unique(t$sg)){
  cor_eff<-c()
  metrics<-c()
  sg<-c()
  g<-t%>% filter(sg==s) 
  for (i in 2:6){
    cor_eff<- c(cor_eff, cor.test(g[, "Fitness"] %>% unlist(), g[, i] %>% unlist(), method="spearman")$estimate)
    metrics<-c(metrics, colnames(g)[i])
    sg<-c(sg, s)
  }
  
  g<-g %>%gather(metrics, values, 2:6)
  g$values[g$metrics %in% c("EvoEF", "arDCA")]<- -1 * g$values[g$metrics %in% c("EvoEF", "arDCA")]
  b<- right_join(
    left_join(g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
                filter(rank_n>=0.95, rank_p>=0.95) %>% group_by(sg, metrics) %>% summarise(n=n()),
              g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
                filter(rank_n>=0.95) %>% group_by(sg, metrics) %>% summarise(t=n())) %>% mutate(overlap_5=n/t) %>%
      select(sg, metrics, overlap_5),
    data.frame(sg, metrics, cor_eff)) %>% mutate(datatype="SAV_SaCas9")
  tab<-bind_rows(tab, b)
}
SAV1<-"THRNTNNQDAAGHRGKNKH"
SAV2<-"THRNTNNQDRQGYAAKNKH"

### look for consensus
t<-t %>% gather(predictor, values, 2:5)
t$values<-as.numeric(t$values)
t$values[t$predictor %in% c("EvoEF", "arDCA")]<- -1 * t$values[t$predictor %in% c("EvoEF", "arDCA")]
### try to show the selections in graph
### plot 1) boxplot of the predicted top 50 fitness and plot 2) no. of diff from WT
t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  group_by(AACombo, sg) %>% summarise(c=n())  

  t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  filter(AACombo %in% c(SAV1, SAV2))


### Cas13d
t<-read_excel(skip=1, "tableS1.xlsx", sheet = 5)
t<-t %>% gather(sg, Fitness, 3)
t<-t [, c( "AACombo", "EvoEF","EVmutation", "Prot_bert-Naive", "MSA-Naive", "arDCA","sg","Fitness")]
s<-unique(t$sg)
cor_eff<-c()
metrics<-c()
sg<-c()
g<-t[!is.na(t$AACombo), ]
for (i in 2:6){
  cor_eff<- c(cor_eff, cor.test(g[, "Fitness"] %>% unlist(), g[, i] %>% unlist(), method="spearman")$estimate)
  metrics<-c(metrics, colnames(g)[i])
  sg<-c(sg, s)
}

g<-g %>%gather(metrics, values, 2:6)
g$values[g$metrics %in% c("EvoEF", "arDCA")]<- -1 * g$values[g$metrics %in% c("EvoEF", "arDCA")]
b<- right_join(
  left_join(g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
              filter(rank_n>=0.95, rank_p>=0.95) %>% group_by(sg, metrics) %>% summarise(n=n()),
            g %>% group_by(sg, metrics) %>% mutate(rank_n=percent_rank(Fitness), rank_p=percent_rank(values)) %>%
              filter(rank_n>=0.95) %>% group_by(sg, metrics) %>% summarise(t=n())) %>% mutate(overlap_5=n/t) %>%
    select(sg, metrics, overlap_5),
  data.frame(sg, metrics, cor_eff)) %>% mutate(datatype="rfx_CAs13d")

tab<-bind_rows(tab, b)
tab[is.na(tab)]<-0

### look at consensus? 
t<-t %>% gather(predictor, values, 2:6)
t$values<-as.numeric(t$values)
t$values[t$predictor %in% c("EvoEF", "arDCA")]<- -1 * t$values[t$predictor %in% c("EvoEF", "arDCA")]
### try to show the selections in graph
### plot 1) boxplot of the predicted top 50 fitness and plot 2) no. of diff from WT
g<-t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  group_by(AACombo, sg) %>% summarise(c=n())

t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  filter(AACombo %in% c(N2V7, N2V8))
g<-g%>% group_by(sg) %>% mutate(rank_c=rank(-c, ties.method = "min")) 
g<-g %>% mutate(WT_d=stringdist(AACombo, WT))

N2V7<-"IIIDALAAYATNAAYAVNNA"
N2V8<-"IIIDILVEYITNVVYVVNNI"

t %>% group_by(predictor, sg) %>%
  mutate(rank_p=percent_rank(values)) %>% filter(rank_p>=0.95) %>%
  filter(AACombo %in% c(N2V7, N2V8))

tab<-read_csv("Zeroshot_performance_summary.csv")
tab$metrics[grepl("esm_msa1", tab$metrics)]<-"MSA"
tab$metrics[grepl("EvoEF", tab$metrics)]<-"EvoEF"
tab<- tab %>% filter(datatype!="SaCas9_888889909")
### plot out how often the predictor rank #1 in each sgRNA
### remove overlap_5 == 0 as it doesn't help
### arDCA sometimes has the wrong sign on cor-eff
#p<-tab[(tab$metrics %in% c("arDCA", "EvoEF") & tab$cor_eff < 0) | (!tab$metrics %in% c("arDCA", "EvoEF") & tab$cor_eff > 0) , ]
# if tab[(tab$metrics %in% c("arDCA", "EvoEF") & tab$cor_eff >0 ] cor_eff = 0
tab$cor_eff[tab$metrics %in% c("arDCA", "EvoEF") & tab$cor_eff >0 ]<-0
tab$cor_eff[(!tab$metrics %in% c("arDCA", "EvoEF") & tab$cor_eff < 0)]<-0
### so now any cor with the wrong sign are purge ###

p1<-tab %>% gather(e, values, 3:4) 
p1$e[p1$e=="cor_eff"]<-"Spearman's rho"
p1$e[p1$e=="overlap_5"]<-"Overlap_5"

p1<-p1%>% ggplot(aes(x=metrics, y=values, colour=datatype, group=metrics)) +
  geom_boxplot(position = position_dodge(width = 0.6))+
  geom_point(position = position_jitterdodge(dodge.width= 0.6, jitter.width = 0.3))+
  geom_hline(data= p1%>% filter(e=="Overlap_5"), aes(yintercept = 0.05), linetype=2)+
  theme_publication+
  facet_grid(.~e, scales = "free_x")+
  coord_flip()+
  scale_color_npg()+
  theme(axis.title = element_blank(),
        legend.position = "right")

tab %>% group_by(metrics) %>% summarise(mean(cor_eff))
tab %>% group_by(metrics) %>% summarise(mean(overlap_5))


max_tab<-bind_rows(
left_join(tab %>% group_by(sg) %>% summarise(overlap_5=max(overlap_5)) %>% filter(overlap_5>0),
          tab) %>% group_by(metrics) %>% summarise(max_o=n()) %>% gather(metric, max, 2),

left_join(tab %>% mutate(cor_eff=abs(cor_eff))  %>% group_by(sg) %>% summarise(cor_eff=max(cor_eff)),
 tab  %>% mutate(cor_eff=abs(cor_eff))) %>% group_by(metrics) %>% 
  summarise(max_c=n())  %>% gather(metric, max, 2))



max_tab<-left_join(expand_grid(metrics= unique(tab$metrics), metric=c("max_o", "max_c")),
          max_tab)
max_tab$max[is.na(max_tab$max)]<-0
colnames(max_tab)<-c("predictor", "metrics", "max")

max_tab$predictor[max_tab$predictor=="EvoEF"]<-"EvoEF"

max_tab$metrics[max_tab$metrics=="max_o"]<-"Overlap_5"
max_tab$metrics[max_tab$metrics=="max_c"]<-"Spearman's rho"


p2<-max_tab %>% ggplot(aes(y=max, x=predictor, fill=predictor))+
  geom_bar(stat = "identity")+
  facet_grid(.~metrics)+
  theme_publication+
  coord_flip()+
  theme(strip.background = element_blank())+
  labs(y="Number of times as best predictor")+
  scale_fill_manual(values=c("MSA"="darkseagreen", "EvMutation"="#4DBBD5FF", 
"prot_bert-Naive"="#00A087FF", "EvoEF"="#3C5488FF", "arDCA"="deeppink"))+
  theme(legend.position = "none",
        axis.title = element_blank())

library(patchwork)
p1 + p2 + plot_layout(ncol = 1)
ggsave("Fig1A.svg", plot=p1 + p2 + plot_layout(ncol = 1), width=8, height=3)

### plot mutation distance
SpCas9<-read_csv("SpCas9_zeroshotPred.csv")
SaCas9<-read_csv("SaCas9_Wed-PI_zeroshotPred.csv")
library(stringdist)
WT<-"RQKETQKR"
SpCas9$mutdis<-stringdist(SpCas9$AACombo, WT)
WT<-"LNANNLLR"
SaCas9$mutdis<-stringdist(SaCas9$AACombo, WT)
colnames(SpCas9)[4:8]<-c("EvoEF", "EVmutation", "prot_bert-Naive" , "MSA", "arDCA")
colnames(SaCas9)[5:9]<-c("EvoEF", "EVmutation", "prot_bert-Naive" , "MSA", "arDCA")

SpCas9<-SpCas9 %>% gather(predictor, predictedfitness, 4:8)
SaCas9<-SaCas9 %>% gather(predictor, predictedfitness, 5:9)

p1<-bind_rows(
SpCas9 %>% group_by(predictor) %>% mutate(rank_p=percent_rank(predictedfitness), datatype="SpCas9") %>% filter(rank_p>=0.95) %>% select(datatype, predictor, mutdis),
SaCas9 %>% group_by(predictor) %>% mutate(rank_p=percent_rank(predictedfitness), datatype="SaCas9_WED-PI") %>% filter(rank_p>=0.95)%>% select(datatype, predictor, mutdis)) 
p1$predictor<-factor(p1$predictor, c("prot_bert-Naive", "EVmutation", "EvoEF",  "MSA", "arDCA"))
p<-p1%>%
  ggplot(aes(x=predictor, y=mutdis))+
  geom_boxplot()+
  facet_grid(.~datatype)+
  coord_flip()+
  theme_publication+
  labs(y="AA distance from WT")
ggsave(p, filename = "figS1AAdis.svg", width=8, height=3.5)
