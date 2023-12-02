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
### input data supplementary table S5
q<-read_xlsx(skip=1, "tableS5.xlsx", sheet=1)
colnames(q)[2]<-"sg"
q$sampling[q$sampling=="EVmutation_CSTopCSTop"]="EVmutation_TopVIP"
q$sampling[q$sampling=="arDCA_CSTopCSTop"]="arDCA_TopVIP"
q$sampling[q$sampling=="EVmutation+arDCA_CS2TopCSTop"]="EVmutation+arDCA_TopVIP+CLADE2.0"

q$sampling<-factor(q$sampling, levels=c("augmentedEVmutation", "augmentedVAE", "seqdesign", "efficient-evolution", "arDCA_TopVIP", "EVmutation_TopVIP", "EVmutation+arDCA_TopVIP+CLADE2.0"))
q$protein<-factor(q$protein, levels=c("SpCas9", "SaCas9_1296", "SaCas9_8000", "evoAPOBEC", "GFP", "PABP", "PhoQ", "UBE"))
q %>% ggplot(aes(x=sampling, y=NDCG, group=sampling)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width =0.05))+
  theme_publication+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylim(0.9, 1)

q %>% ggplot(aes(x=sampling, y=`top 50 identified`, group=sampling)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width =0.05))+
  theme_publication+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

q %>% ggplot(aes(x=sampling, y=selectivity, group=sampling)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width =0.05))+
  theme_publication+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

q %>% ggplot(aes(x=sampling, y=t5, group=sampling)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width =0.05))+
  theme_publication+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

t<-q %>% gather(metric, values, 4:10)
#t<-t%>% filter(!grepl("CLADE2EV", sampling))

t<- t %>% filter(metric!= "top 50 in training data")
t$metric<-factor(t$metric, levels=c("top 50 identified", "top 1% identified", "enrichment", "selectivity", "NDCG", "Spearnman's rho"))

p1<-q 

a1<-p1%>% 
  ggplot(aes(x=selectivity, y=`top 1% identified`/(`#total variants`/100), colour=sampling)) +
  geom_point(size=2)+
  facet_wrap(~protein, scales = "free", ncol=4)+
  labs(x="selectivity", y="discovery rate of top 1% variants")+ theme_publication+
  scale_color_aaas()
ggsave(plot=a1, filename = "benchmark_top1vsselectivity.svg", width=8, height=5)


### ranking
r<-q %>% group_by(sampling)%>% summarise(cases=n())

r<-left_join(r, q %>% group_by(protein, sg) %>%
  mutate(rank_top50=rank(-`top 50 identified`, ties.method = "min")) %>%
  filter(rank_top50==1) %>% group_by(sampling)%>% summarise(top_50=n()) %>% arrange(desc(top_50))
)
r<-left_join(r, q  %>% group_by(protein, sg) %>%
  mutate(rank_top50=rank(-selectivity, ties.method = "min")) %>%
  filter(rank_top50==1) %>% group_by(sampling)%>% summarise(selectivity=n()) %>% arrange(desc(selectivity)) 
)
r<-left_join(r,q  %>% group_by(protein, sg) %>%
  mutate(rank_top50=rank(-NDCG, ties.method = "min")) %>%
  filter(rank_top50==1) %>% group_by(sampling)%>% summarise(NDCG=n()) %>% arrange(desc(NDCG)) 
)
r<-left_join(r,q %>% group_by(protein, sg) %>%
  mutate(rank_top50=rank(-enrichment, ties.method = "min")) %>%
  filter(rank_top50==1) %>% group_by(sampling)%>% summarise(enrichment=n()) %>% arrange(desc(enrichment)) 
)
r<-left_join(r,q %>% group_by(protein, sg) %>%
               mutate(rank_top1p=rank(-`top 1% identified`, ties.method = "min")) %>%
               filter(rank_top1p==1) %>% group_by(sampling)%>% summarise(top_1_percent=n()) %>% arrange(desc(top_1_percent)) 
)

r[is.na(r)]<-0
r<-r %>% gather(metric, freq, 3:7) 
r$metric<- factor(r$metric, levels=c("top_50", "top_1_percent",  "enrichment", "selectivity", "NDCG"))
a1<-r %>% 
  ggplot(aes(x=sampling, y=freq))+
  geom_bar(aes(fill=sampling), stat="identity")+
  facet_grid(.~metric)+
  coord_flip()+
  labs(x="No. of times as best predictor")+
  theme_publication+
  scale_fill_aaas()+
  scale_y_continuous(breaks= seq(0, 10, 2))
ggsave(plot=a1, filename = "benchmark_best_compare.svg", width=8.5, height=3)


