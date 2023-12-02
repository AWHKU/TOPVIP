library(ggplot2)
library(tidyverse)
library(ggpmisc)
library(ggseqlogo)
library(readxl)
library(Cairo)
library(stringdist)
library(ggsci)
library(ggrepel)
library(tidyverse)
library(readxl)
library(patchwork)
library(ggnewscale)
theme_publication<-theme_classic()+
  theme(text = element_text(size=12, family = "Arial"),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        panel.background = element_rect(colour="black", linewidth=0.5),
        legend.position = "bottom",
        strip.background = element_blank())

### compared with previous 20% random data
setwd("/media/achu/新增磁碟區/CLADE/2023_12_01_fdraft/data")
w_o<-read_csv("/media/achu/新增磁碟區/CLADE/2023_12_01_fdraft/data/performance_20_percent.csv")
#### 1-round comparison
w<-read_csv("/media/achu/新增磁碟區/CLADE/2023_12_01_fdraft/data/1-round_MLDE_performance.csv")
### divide into datasets ###
a1<-w %>% filter(datatype=="SpCas9") %>% filter(no<96, sampling_r1!= "random+WT") %>%
  ggplot(aes(x=sampling_r1, y=performance, color=predictor, fill=sampling, group=interaction(no, sampling_r1, predictor)))+
  geom_hline(data=w_o %>% filter(datatype =="SpCas9", !metrics %in% c("top 20", "top 50")), aes(yintercept=performance), linetype=2, colour="purple")+
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.05))+
  geom_boxplot(position = position_dodge(width=0.8))+
  facet_grid(metrics~no, scales = "free")+
  theme_publication+
  scale_color_manual(values=c("Evmutation"="#4DBBD5FF", "EvoEF"="#3C5488FF", 
                              "MSA"="darkseagreen","arDCA"="deeppink", "random"="black"))+
  scale_fill_manual(values=c("CS"="grey50", "direct-ranking"="white", "random"="grey20"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        legend.box="vertical", 
        legend.box.margin = margin(0,0,0,0),
        strip.placement = "outside")+
  labs(title="SpCas9")

a2<-w %>% filter(datatype=="SaCas9_WED-PI") %>% filter(no<96, sampling_r1!= "random+WT") %>%
  ggplot(aes(x=sampling_r1, y=performance, color=predictor, fill=sampling, group=interaction(no, sampling_r1, predictor)))+
  geom_hline(data=w_o %>% filter(datatype =="SaCas9_WED-PI", !metrics %in% c("top 20", "top 50")), aes(yintercept=performance), linetype=2, colour="purple")+
  geom_boxplot(position = position_dodge(width=0.8))+
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.05))+
  facet_grid(metrics~no, scales = "free")+
  theme_publication+
  scale_color_manual(values=c("Evmutation"="#4DBBD5FF", "EvoEF"="#3C5488FF", 
                              "MSA"="darkseagreen", "arDCA"="deeppink", "random"="black"))+
  scale_fill_manual(values=c("CS"="grey50", "direct-ranking"="white", "random"="grey20"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        legend.box="vertical", 
        legend.box.margin = margin(0,0,0,0),
        strip.placement = "outside")+
  labs(title="SaCas9_WED-PI")

library(patchwork)
a1/a2
#ggsave(plot=a1, filename = "Fig_s1.svg", width=8, height=10)
#ggsave(plot=a2, filename = "Fig_s2.svg", width=8, height=10)

### compared to random_NO_NO, how much have each strategy improve the metrics?
improve_tab<- left_join(w %>% filter(predictor!="random"),
w %>% filter(sampling_strategy=="random_NO_NO") %>% mutate(random_NO_NO=performance) %>%
  select(no, datatype, sg, metrics, random_NO_NO)) %>% mutate(improvement=(performance-random_NO_NO)) %>%
  group_by(datatype, predictor, sampling_strategy, metrics) %>% mutate(m_improve=mean(improvement, na.rm=T))
improve_tab$sampling_strategy<-factor(improve_tab$sampling_strategy, levels=c("direct-ranking_NO_NO", "CLADE_NO_NO", "CLADE_YES_NO","CLADE_NO_YES", "CLADE_YES_YES", "random_NO_NO", "random_NO_YES"  ))

improve_tab$sampling_strategy_rename<-gsub("CLADE", "CS", improve_tab$sampling_strategy)
improve_tab %>% ggplot(aes(x=interaction(sampling_strategy_rename, datatype), y=m_improve, fill=sampling_strategy))+
  geom_bar(stat="unique", position="identity")+
 geom_point(aes(y=improvement), size=0.7, colour="grey50")+
  facet_grid(metrics~predictor, scales = "free")+
  theme_publication+
  labs(y="Improvement", x="sampling_strategy")+
  theme(strip.placement = "outer")+
  theme(axis.text.x = element_text(angle = 90))
### t.test
q1<-left_join(
  left_join(
    left_join(
      left_join(
        left_join(
          w %>% filter(sampling_strategy=="random_NO_NO") %>% mutate(random_NO_NO=performance) %>% select(no, datatype, sg, metrics, random_NO_NO) , 
          w %>% filter(sampling_strategy=="direct-ranking_NO_NO") %>% mutate(direct_ranking=performance) %>% select(no, datatype, predictor,  sg, metrics,direct_ranking)),
        w %>% filter(sampling_strategy=="CLADE_NO_NO") %>% mutate(CLADE_NO_NO=performance) %>% select(no, datatype, predictor,  sg, metrics, CLADE_NO_NO) ),
      w %>% filter(sampling_strategy=="CLADE_NO_YES") %>% mutate(CLADE_NO_YES=performance) %>% select(no, datatype, predictor,  sg, metrics, CLADE_NO_YES) ),
    w %>% filter(sampling_strategy=="CLADE_YES_NO") %>% mutate(CLADE_YES_NO=performance) %>% select(no, datatype, predictor,  sg, metrics, CLADE_YES_NO) ),  
  w %>% filter(sampling_strategy=="CLADE_YES_YES") %>% mutate(CLADE_YES_YES=performance) %>% select(no, datatype, predictor,  sg, metrics, CLADE_YES_YES)  )

q1<-q1 %>% filter(no<96)

s1<-q1 %>% select(no, metrics, direct_ranking, CLADE_NO_NO, CLADE_YES_NO, CLADE_NO_YES, CLADE_YES_YES)
colnames(s1)[3:7]<-c("direct-ranking", "CLADE","CLADE+AS","CLADE+WT","CLADE+AS+WT")
s1<-s1%>%
  gather(sampling_strategy, improvement, 3:7)
s1$sampling_strategy<-factor(s1$sampling_strategy, levels=c("direct-ranking", "CLADE","CLADE+AS","CLADE+WT","CLADE+AS+WT","random", "random+WT"))

s1$sampling_strategy_rename<-gsub("CLADE", "CS", s1$sampling_strategy)
s1$sampling_strategy_rename<-factor(s1$sampling_strategy_rename, levels=c("direct-ranking", "CS","CS+AS","CS+WT","CS+AS+WT","random", "random+WT"))

b1<-s1 %>% filter(metrics=="NDCG") %>%
  ggplot(aes(x=sampling_strategy_rename, y=improvement, colour=sampling_strategy))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(group=no), position = position_dodge(width=0.7), size=0.5)+
  theme_publication+
  theme(strip.placement = "outer")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")+
  labs(y="NDCG")

b2<-s1 %>% filter(metrics=="enrichment") %>%
  ggplot(aes(x=sampling_strategy_rename, y=improvement, colour=sampling_strategy))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(group=no), position = position_dodge(width=0.7), size=0.5)+
  theme_publication+
  theme(strip.placement = "outer")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")+
  labs(y="enrichment")

b3<-s1 %>% filter(metrics=="rho") %>%
  ggplot(aes(x=sampling_strategy_rename, y=improvement, colour=sampling_strategy))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(group=no), position = position_dodge(width=0.7), size=0.5)+
  theme_publication+
  theme(strip.placement = "outer")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")+
  labs(y="Spearman's rho")


b4<-s1 %>% filter(metrics=="selectivity") %>%
  ggplot(aes(x=sampling_strategy_rename, y=improvement, colour=sampling_strategy))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(group=no), position = position_dodge(width=0.7), size=0.5)+
  theme_publication+
  theme(strip.placement = "outer")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")+
  labs(y="selectivity")

library(patchwork)
ggsave((b1+b2)/(b3+b4), file="improvement_t.test.svg", width=6, height=6)

q1<-left_join(
left_join(
left_join(
left_join(
left_join(
w %>% filter(sampling_strategy=="random_NO_NO") %>% mutate(random_NO_NO=performance) %>% select(no, datatype, sg, metrics, random_NO_NO) , 
w %>% filter(sampling_strategy=="direct-ranking_NO_NO") %>% mutate(direct_ranking=performance) %>% select(no, datatype, predictor,  sg, metrics,direct_ranking)),
w %>% filter(sampling_strategy=="CLADE_NO_NO") %>% mutate(CLADE_NO_NO=performance) %>% select(no, datatype, predictor,  sg, metrics, CLADE_NO_NO) ),
w %>% filter(sampling_strategy=="CLADE_NO_YES") %>% mutate(CLADE_NO_YES=performance) %>% select(no, datatype, predictor,  sg, metrics, CLADE_NO_YES) ),
w %>% filter(sampling_strategy=="CLADE_YES_NO") %>% mutate(CLADE_YES_NO=performance) %>% select(no, datatype, predictor,  sg, metrics, CLADE_YES_NO) ),  
w %>% filter(sampling_strategy=="CLADE_YES_YES") %>% mutate(CLADE_YES_YES=performance) %>% select(no, datatype, predictor,  sg, metrics, CLADE_YES_YES)  )

q1<-q1 %>% filter(no<96)
### NDCG - among datatype
improve_tab<- left_join(w %>% filter(predictor!="random"),
                        w %>% filter(sampling_strategy=="random_NO_NO") %>% mutate(random_NO_NO=performance) %>%
                          select(no, datatype, sg, metrics, random_NO_NO)) %>% mutate(improvement=(performance-random_NO_NO)) %>%
  group_by( predictor, sampling_strategy, metrics, no) %>% mutate(m_improve=mean(improvement, na.rm=T))
improve_tab$sampling_strategy<-factor(improve_tab$sampling_strategy, levels=c("direct-ranking_NO_NO", "CLADE_NO_NO", "CLADE_YES_NO","CLADE_NO_YES", "CLADE_YES_YES", "random_NO_NO", "random_NO_YES"  ))
improve_tab<-left_join(improve_tab, w_o %>% mutate(p20=performance) %>% select(datatype, metrics, p20))
improve_tab<-improve_tab %>% mutate(ub=p20-random_NO_NO)
improve_tab$sampling_strategy_rename<-gsub("CLADE", "CS", improve_tab$sampling_strategy)

p1<-improve_tab%>% filter(metrics == "NDCG", no <=48) %>% ggplot(aes(x=sampling_strategy, y=m_improve, fill=sampling_strategy))+
  geom_hline(yintercept = 0)+
  geom_hline(data=improve_tab %>% filter(metrics == "NDCG", no <=48) %>% group_by(no, predictor) %>% summarise(ub=mean(ub)), aes(yintercept=ub), color="purple", linetype=2 )+
  geom_bar(stat="identity", position="identity")+
  geom_point(aes(y=improvement), size=0.7, colour="grey50")+
  facet_grid(predictor~no)+
  theme_publication+
  labs(y="Improvement on NDCG", x="sampling_strategy")+
  scale_x_discrete(labels=c("direct-ranking", "CS", "CS+AS", "CS+WT", "CS+AS+WT"))+
  theme(strip.placement = "outer")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
legend.position = "none")
ggsave(p1, filename = "/media/achu/新增磁碟區/CLADE/Figures/improvement_NDCG.svg", height=5, width=4)

p1<-improve_tab%>% filter(metrics == "enrichment", no <=48) %>% ggplot(aes(x=sampling_strategy, y=m_improve, fill=sampling_strategy))+
  geom_hline(yintercept = 0)+
  geom_hline(data=improve_tab %>% filter(metrics == "enrichment", no <=48) %>% group_by(no, predictor) %>% summarise(ub=mean(ub)), aes(yintercept=ub), color="purple", linetype=2 )+
  geom_bar(stat="identity", position="identity")+
  geom_point(aes(y=improvement), size=0.7, colour="grey50")+
  facet_grid(predictor~no)+
  theme_publication+
  labs(y="Improvement on enrichment", x="sampling_strategy")+
  scale_x_discrete(labels=c("direct-ranking", "CS", "CS+AS", "CS+WT", "CS+AS+WT"))+
  theme(strip.placement = "outer")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")
ggsave(p1, filename = "/media/achu/新增磁碟區/CLADE/Figures/improvement_enrichment.svg", height=5, width=3.8)

p1<-improve_tab%>% filter(metrics == "selectivity", no <=48) %>% ggplot(aes(x=sampling_strategy, y=m_improve, fill=sampling_strategy))+
  geom_hline(yintercept = 0)+
  geom_bar(stat="identity", position="identity")+
  geom_hline(data=improve_tab %>% filter(metrics == "selectivity", no <=48) %>% group_by(no, predictor) %>% summarise(ub=mean(ub)), aes(yintercept=ub), color="purple", linetype=2 )+
  geom_point(aes(y=improvement), size=0.7, colour="grey50")+
  facet_grid(predictor~no, scales = "free")+
  theme_publication+
  labs(y="Improvement on selectivity", x="sampling_strategy")+
  scale_x_discrete(labels=c("direct-ranking", "CS", "CS+AS", "CS+WT", "CS+AS+WT"))+
  theme(strip.placement = "outer")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")
ggsave(p1, filename = "/media/achu/新增磁碟區/CLADE/Figures/improvement_selectivity.svg", height=5, width=3.8)

### 1-2-4 round sampling 
### boxplot 1,2,4 sampling ###
test_p1<-read_csv("/media/achu/新增磁碟區/CLADE/2023_12_01_fdraft/data/1-2-4-round_sampling_performance.csv")
test_p1$no<- factor(test_p1$no, levels=c("48", "24+24", "12x4"))
test_p1$sg<-factor(test_p1$sg, levels = c("Sg5", "Sg8", 
                              "ACAAGT", "GGTGGT", "TGGAGT"))
test_p1$predictor<-factor(test_p1$predictor, c("Evmutation", "EvoEF", "MSA",  "arDCA", "random"))
test_p1$sampling_r2_rename<-gsub("CLADE", "CS", test_p1$sampling_r2)

a1<-test_p1 %>% filter(datatype=="SpCas9") %>%
  ggplot(aes(x=interaction(sampling_r1, sampling_r2_rename), y=performance, fill=predictor, color=predictor, 
             group=interaction(predictor, sampling_r1, sampling_r2)))+
  geom_hline(data=w_o %>% filter(datatype=="SpCas9", metrics != "top 20"), aes(yintercept=performance), linetype=2, colour="purple")+
  geom_boxplot(position = position_dodge(width=0.8), outlier.colour = NA)+
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5), size=0.7, colour="grey70")+
  facet_grid(metrics~no, scales = "free", space = "free_x")+
  scale_fill_manual(values=c("Evmutation"="#4DBBD5FF", "EvoEF"="#3C5488FF", 
                                "MSA"="darkseagreen", "arDCA"="deeppink", "random"="black"))+
  scale_color_manual(values=c("Evmutation"="#4DBBD5FF", "EvoEF"="#3C5488FF", 
                           "MSA"="darkseagreen", "arDCA"="deeppink", "random"="black")) +
  theme_publication+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="sampling strategy", title="SpCas9")

a2<-test_p1 %>% filter(datatype=="SaCas9_WED-PI") %>%
  ggplot(aes(x=interaction(sampling_r1, sampling_r2_rename), y=performance, fill=predictor, color=predictor, 
             group=interaction(predictor, sampling_r1, sampling_r2)))+
  geom_hline(data=w_o %>% filter(datatype=="SaCas9_WED-PI", metrics != "top 20"), aes(yintercept=performance), linetype=2, colour="purple")+
  geom_boxplot(position = position_dodge(width=0.8), outlier.colour = NA)+
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5), size=0.7, colour="grey70")+
  facet_grid(metrics~no, scales = "free", space = "free_x")+
  scale_fill_manual(values=c("Evmutation"="#4DBBD5FF", "EvoEF"="#3C5488FF", 
                             "MSA"="darkseagreen", "random"="black", "arDCA"="deeppink"))+
  scale_color_manual(values=c("Evmutation"="#4DBBD5FF", "EvoEF"="#3C5488FF", 
                              "MSA"="darkseagreen", "random"="black", "arDCA"="deeppink")) +
  theme_publication+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="sampling strategy", title="SaCas9_WED-PI")

#ggsave(plot=a1, filename = "FigS3.svg", width=8, height=10)
#ggsave(plot=a2, filename = "FigS4.svg", width=8, height=10)


### scatter plot on NDCG + selectivity?
test_t1<-read_csv("/media/achu/新增磁碟區/CLADE/2023_12_01_fdraft/data/1-2-4-round_scatter_plot.csv")
test_t1$predictor<-factor(test_t1$predictor, c("EVmutation", "EvoEF", "MSA",  "arDCA"))

a1<-test_t1%>%
  ggplot(aes(x=selectivity, y=NDCG, colour=no, shape=predictor))+
  geom_vline(xintercept = 0.05, color="purple", linetype=2)+
  geom_hline(yintercept = 0.950, color="purple", linetype=2)+
  geom_point(size=2)+
  geom_label_repel(data=t1 %>% filter(NDCG>0.945, selectivity>0.3), aes(label=paste0(predictor, "_", sampling_r1)), size=2.5)+
  ylim(0.87, 0.97)+
  xlim(0.0, 0.4)+
  theme_publication+
  theme(legend.box = "vertical")

a2<-test_t1%>%
  ggplot(aes(x=selectivity, y=t50, colour=no, shape=predictor))+
  geom_vline(xintercept = 0.05, color="purple", linetype=2)+
  geom_hline(yintercept = 21.3, color="purple", linetype=2)+
  geom_point(size=2)+
  geom_label_repel(data=t1 %>% filter(t50>23, selectivity>0.32), aes(label=paste0(predictor, "_", sampling_r1)), size=2.5, nudge_y = 1.5)+
  ylim(0, 26)+
  xlim(0.0, 0.4)+
  theme_publication+
  labs(y="Top 50 identified")+
  theme(legend.box = "vertical")

library(patchwork)
ggsave(plot=a1+a2, filename = "predictor_samplingr1_NDCG.svg", width=7, height=4.5)
# base on the best predictors, what kind of sampling scheme is the best?
test_t<-read_csv("/media/achu/新增磁碟區/CLADE/2023_12_01_fdraft/data/1-2-4-round_scatter_plot_fulltab.csv")
test_t$sg<-factor(test_t$sg, levels = c("Sg5", "Sg8", 
                                                "ACAAGT", "GGTGGT", "TGGAGT"))
test_t$predictor<-factor(test_t$predictor, c("EVmutation", "EvoEF", "MSA",  "arDCA", "random"))

t1<-test_t %>% filter(predictor=="arDCA", round==4) %>%  group_by(datatype, no, predictor, sampling_r1, sampling_r2, round, total_n) %>%
  summarise(NDCG=mean(NDCG), selectivity=mean(selectivity), t50=mean(`top 50`))
t1$sampling_r1[t1$sampling_r1=="YES_YES"]<-"AS+WT"
t1$sampling_r1[t1$sampling_r1=="NO_YES"]<-"+WT"
t1$sampling_r1[t1$sampling_r1=="NO_NO"]<-""
t1$sampling_r1[t1$sampling_r1=="YES_NO"]<-"AS"
t1$sampling_r2_rename<-gsub("CLADE", "CS", t1$sampling_r2)
a1<-t1%>%
  ggplot(aes(x=selectivity, y=t50, shape=predictor))+
  geom_point(size=2, shape=3)+
  geom_label_repel(data=t1 %>% filter(t50>24, selectivity>0.4), aes(label=paste0(sampling_r1, "_", sampling_r2_rename)), size=2.5, min.segment.length = 0.3, nudge_y = 1)+
  ylim(10, 30)+
  xlim(0.1, 0.5)+
  theme_publication+
  labs(y="Top 50 identified", title = "arDCA")

t1<-test_t %>% filter(predictor=="EVmutation", round==4) %>%  group_by(datatype, no, predictor, sampling_r1, sampling_r2, round, total_n) %>%
  summarise(NDCG=mean(NDCG), selectivity=mean(selectivity), t50=mean(`top 50`))
t1$sampling_r1[t1$sampling_r1=="YES_YES"]<-"AS+WT"
t1$sampling_r1[t1$sampling_r1=="NO_YES"]<-"+WT"
t1$sampling_r1[t1$sampling_r1=="NO_NO"]<-""
t1$sampling_r1[t1$sampling_r1=="YES_NO"]<-"AS"

t1$sampling_r2_rename<-gsub("CLADE", "CS", t1$sampling_r2)
a2<-t1%>%
  ggplot(aes(x=selectivity, y=t50, shape=predictor, color=datatype))+
  geom_point(size=2, shape=16)+
  geom_label_repel(data=t1 %>% filter(t50>23, selectivity>0.35), aes(label=paste0(sampling_r1, "_", sampling_r2_rename)), size=2.5, min.segment.length = 0.1)+
  ylim(10, 30)+
  xlim(0.1, 0.5)+
  theme_publication+
  labs(y="Top 50 identified", title = "EVmutation")
ggsave(plot=a1+a2, filename = "predictor_samplingr2_compare.svg", width=7, height=3.5)

### Compare
### CLADETopCLADETop or CLADETopTopTop for arDCA and EVmutation
###58libD21 + 56libD18
### show sgRNA empirical ###
Sa_Fit_8000<- read_csv("SaCas9_8888889909_Fitness.csv")
a1<-Sa_Fit_8000 %>% 
  ggplot(aes(x=DTp56libD18, y=DTp58libD21))+
  geom_point(alpha=0.05)+
  geom_abline(slope = 1)+
  geom_point(data=Sa_Fit_8000 %>% mutate(rank_1=rank(-DTp56libD18, ties.method = "min")) %>% filter(rank_1<=50), color="blue", alpha=0.3)+
  geom_point(data=Sa_Fit_8000 %>% mutate(rank_2=rank(-DTp58libD21, ties.method = "min")) %>% filter(rank_2<=50), color="red", alpha=0.3)+
  geom_point(data=Sa_Fit_8000 %>% filter(AACombo=="NAL"), color="orange")+
  geom_text_repel(data=Sa_Fit_8000 %>% mutate(rank_1=rank(-DTp56libD18, ties.method = "min"), p1=substring(AACombo, 1, 1)) %>% filter(rank_1<=50, p1=="G"), aes(label=AACombo), max.overlaps = 100, size=3)+
  geom_text_repel(data=Sa_Fit_8000 %>% mutate(rank_2=rank(-DTp58libD21, ties.method = "min"), p1=substring(AACombo, 1, 1)) %>% filter(rank_2<=50, p1=="G"), aes(label=AACombo), max.overlaps = 100, size=3)+
  theme_publication+
  xlim(1, 3)+
  ylim(1, 6)+
  labs(title = "Empirical")
ggsave(plot=a1, filename = "Fig_3c.svg", width=3, height=3)

### find training datas ###

SaCas9_888889909_prediction<- read_csv("data/SaCas9_888889909_prediction.csv")
v_train<- read_csv("data/SaCas9_888889909_traiing_data.csv")
t<-SaCas9_888889909_prediction %>%   mutate(samplename=paste0(sg, "_", predictor, "_", sampling_r2) ) %>% 
  select(AACombo, sg, `InTrainingData?`, PredictedFitness, Fitness_n, rank_n, samplename, predictor, sampling_r1, sampling_r2) 

f<-left_join(left_join(left_join(t %>% group_by(samplename) %>% mutate(prank_p=percent_rank(PredictedFitness), prank_n=percent_rank(Fitness_n)) %>% filter(prank_n>=0.95 & `InTrainingData?`=="YES")
                                 %>% group_by(samplename) %>% summarise(n=n()), 
                                 t %>% filter(`InTrainingData?`=="YES") %>% group_by(samplename) %>% summarise(s=n())) %>% mutate(selecitivity=n/s),
                       t %>% group_by(samplename) %>% mutate(rank_p=rank(-PredictedFitness, ties.method  ="min")) %>%
                         filter(rank_n<=50, rank_p<=50) %>%
                         group_by(samplename) %>% summarise(top_50=n())),
             t %>% filter(rank_n<=50, `InTrainingData?`=="YES") %>%
               group_by(samplename) %>% summarise(top_50_intrainingdata=n()))

f$sg<-sapply(f$samplename, function(s) strsplit(s, "_")[[1]][1])
f$predictor<-sapply(f$samplename, function(s) strsplit(s, "_")[[1]][2])
f$sampling<-sapply(f$samplename, function(s) strsplit(s, "_")[[1]][3])
t1<-f
t1$samplename<-gsub("CLADE", "CS", t1$samplename)                
#write_csv(t1, "table1.csv")


#1) arDCA, CLADETopCLADETop
t<-SaCas9_888889909_prediction %>%   mutate(samplename=paste0(sg, "_", predictor, "_", sampling_r2) ) %>% 
  select(AACombo, sg, `InTrainingData?`, PredictedFitness, Fitness_n, rank_n, samplename, predictor, sampling_r1, sampling_r2) 
left_join(left_join(t %>% group_by(samplename) %>% mutate(prank_p=percent_rank(PredictedFitness), prank_n=percent_rank(Fitness_n)) %>% filter(prank_n>=0.95 & `InTrainingData?`=="YES")
          %>% group_by(samplename) %>% summarise(n=n()), 
          t %>% filter(`InTrainingData?`=="YES") %>% group_by(samplename) %>% summarise(s=n())) %>% mutate(n/s),
t %>% group_by(samplename) %>% mutate(rank_p=rank(-PredictedFitness, ties.method  ="min")) %>%
  filter(rank_n<=50, rank_p<=50) %>%
  group_by(samplename) %>% summarise(t=n()))
t %>% filter(rank_n<=50, `InTrainingData?`=="YES") %>%
  group_by(samplename) %>% summarise(v=n())
### only ploting predicted fitness here
t<-t %>% filter(predictor=="arDCA", sampling_r2=="CLADETopCLADETop") %>%
  select(AACombo, PredictedFitness, predictor, sg,  sampling_r1, sampling_r2) %>%
  spread(sg, PredictedFitness)
t<-t%>% mutate(rank1=rank(-DTp56libD18, ties.method  ="min"),
              rank2=rank(-DTp58libD21, ties.method  ="min")) 
arDCA_Top_CTCT<- t %>%   filter(rank1<=50 & rank2<=50) %>% select(AACombo) %>% unlist()
a1<-t %>% ggplot(aes(x=DTp56libD18, y=DTp58libD21, label=AACombo))+
  geom_point(alpha=0.05)+
  geom_point(data=t %>% filter(rank1<=50), colour="cyan", alpha=0.5)+
  geom_point(data=t %>% filter(rank2<=50), colour="deeppink", alpha=0.4)+
  geom_point(data=t %>% filter(AACombo=="NAL"), color="green")+
  geom_point(data= t %>% filter(AACombo %in% v_train$AACombo[v_train$sampling_r2=="CLADETopCLADETop" & v_train$sg=="DTp56libD18" & v_train$predictor=="arDCA"]), pch=1, colour="blue", size=2.5)+
  geom_point(data= t %>% filter(AACombo %in% v_train$AACombo[v_train$sampling_r2=="CLADETopCLADETop" & v_train$sg=="DTp58libD21" & v_train$predictor=="arDCA"]), pch=1, colour="red", size=1.5)+
  geom_text_repel(data= t %>% 
                    filter(rank1<=50 & rank2<=50), 
                  aes(label=AACombo), max.overlaps = 100, size=3, min.segment.length = 0.01)+
  geom_text_repel(data=t %>% filter(AACombo=="NAL"), color="green", size=3.5)+
  xlim(0.5, 1)+
  ylim(0.5, 1.1)+
  theme_publication+
  labs(title="arDCA-+WT_CSTopCSTop", x="DTp56libD18-predictedFitness", y="DTp58libD21-predictedFitness")
ggsave(plot=a1, filename = "Sa888889909_arDCaCTCT.svg", width=3, height=3)

#2) arDCA, CLADETopTopTop
t<-SaCas9_888889909_prediction %>%   mutate(samplename=paste0(sg, "_", predictor, "_", sampling_r2) ) %>% 
  select(AACombo, sg, `InTrainingData?`, PredictedFitness, Fitness_n, rank_n, samplename, predictor, sampling_r1, sampling_r2) 
left_join(left_join(t %>% group_by(samplename) %>% mutate(prank_p=percent_rank(PredictedFitness), prank_n=percent_rank(Fitness_n)) %>% filter(prank_n>=0.95 & `InTrainingData?`=="YES")
                    %>% group_by(samplename) %>% summarise(n=n()), 
                    t %>% filter(`InTrainingData?`=="YES") %>% group_by(samplename) %>% summarise(s=n())) %>% mutate(n/s),
          t %>% group_by(samplename) %>% mutate(rank_p=rank(-PredictedFitness, ties.method  ="min")) %>%
            filter(rank_n<=50, rank_p<=50) %>%
            group_by(samplename) %>% summarise(t=n()))
t %>% filter(rank_n<=50, `InTrainingData?`=="YES") %>%
  group_by(samplename) %>% summarise(v=n())
### only ploting predicted fitness here
t<-t %>% filter(predictor=="arDCA", sampling_r2=="CLADETopTopTop") %>%
  select(AACombo, PredictedFitness, predictor, sg,  sampling_r1, sampling_r2) %>%
  spread(sg, PredictedFitness)
t<-t%>% mutate(rank1=rank(-DTp56libD18, ties.method  ="min"),
               rank2=rank(-DTp58libD21, ties.method  ="min")) 
arDCA_Top_CTTT<- t %>%   filter(rank1<=50 & rank2<=50) %>% select(AACombo) %>% unlist()
a1<-t %>% ggplot(aes(x=DTp56libD18, y=DTp58libD21, label=AACombo))+
  geom_point(alpha=0.05)+
  geom_point(data=t %>% filter(rank1<=50), colour="cyan", alpha=0.5)+
  geom_point(data=t %>% filter(rank2<=50), colour="deeppink", alpha=0.4)+
  geom_point(data=t %>% filter(AACombo=="NAL"), color="green")+
  geom_point(data= t %>% filter(AACombo %in% v_train$AACombo[v_train$sampling_r2=="CLADETopTopTop" & v_train$sg=="DTp56libD18" & v_train$predictor=="arDCA"]), pch=1, colour="blue", size=2.5)+
  geom_point(data= t %>% filter(AACombo %in% v_train$AACombo[v_train$sampling_r2=="CLADETopTopTop" & v_train$sg=="DTp58libD21" & v_train$predictor=="arDCA"]), pch=1, colour="red", size=1.5)+
  geom_text_repel(data= t %>% 
                    filter(rank1<=50 & rank2<=50), 
                  aes(label=AACombo), max.overlaps = 100, size=3, min.segment.length = 0.01)+
  geom_text_repel(data=t %>% filter(AACombo=="NAL"), color="green", size=3.5)+
  xlim(0.5, 1)+
  ylim(0.5, 1.1)+
  theme_publication+
  labs(title="arDCA-+WT_CSTopTopTop", x="DTp56libD18-predictedFitness", y="DTp58libD21-predictedFitness")
ggsave(plot=a1, filename = "Sa888889909_arDCaCTTT.svg", width=3, height=3)

#1) EVmutation, CLADETopCLADETop
t<-SaCas9_888889909_prediction %>%   mutate(samplename=paste0(sg, "_", predictor, "_", sampling_r2) ) %>% 
  select(AACombo, sg, `InTrainingData?`, PredictedFitness, Fitness_n, rank_n, samplename, predictor, sampling_r1, sampling_r2) 
left_join(left_join(t %>% group_by(samplename) %>% mutate(prank_p=percent_rank(PredictedFitness), prank_n=percent_rank(Fitness_n)) %>% filter(prank_n>=0.95 & `InTrainingData?`=="YES")
                    %>% group_by(samplename) %>% summarise(n=n()), 
                    t %>% filter(`InTrainingData?`=="YES") %>% group_by(samplename) %>% summarise(s=n())) %>% mutate(n/s),
          t %>% group_by(samplename) %>% mutate(rank_p=rank(-PredictedFitness, ties.method  ="min")) %>%
            filter(rank_n<=50, rank_p<=50) %>%
            group_by(samplename) %>% summarise(t=n()))
t %>% filter(rank_n<=50, `InTrainingData?`=="YES") %>%
  group_by(samplename) %>% summarise(v=n())
### only ploting predicted fitness here
t<-t %>% filter(predictor=="EVmutation", sampling_r2=="CLADETopCLADETop") %>%
  select(AACombo, PredictedFitness, predictor, sg,  sampling_r1, sampling_r2) %>%
  spread(sg, PredictedFitness)
t<-t%>% mutate(rank1=rank(-DTp56libD18, ties.method  ="min"),
               rank2=rank(-DTp58libD21, ties.method  ="min")) 
EVmutation_Top_CTCT<- t %>%   filter(rank1<=50 & rank2<=50) %>% select(AACombo) %>% unlist()
a1<-t %>% ggplot(aes(x=DTp56libD18, y=DTp58libD21, label=AACombo))+
  geom_point(alpha=0.05)+
  geom_point(data=t %>% filter(rank1<=50), colour="cyan", alpha=0.5)+
  geom_point(data=t %>% filter(rank2<=50), colour="deeppink", alpha=0.4)+
  geom_point(data=t %>% filter(AACombo=="NAL"), color="green")+
  geom_point(data= t %>% filter(AACombo %in% v_train$AACombo[v_train$sampling_r2=="CLADETopCLADETop" & v_train$sg=="DTp56libD18" & v_train$predictor=="EVmutation"]), pch=1, colour="blue", size=2.5)+
  geom_point(data= t %>% filter(AACombo %in% v_train$AACombo[v_train$sampling_r2=="CLADETopCLADETop" & v_train$sg=="DTp58libD21" & v_train$predictor=="EVmutation"]), pch=1, colour="red", size=1.5)+
  geom_text_repel(data= t %>% 
                    filter(rank1<=50 & rank2<=50), 
                  aes(label=AACombo), max.overlaps = 100, size=3, min.segment.length = 0.01)+
  geom_text_repel(data=t %>% filter(AACombo=="NAL"), color="green", size=3.5)+
  xlim(0.5, 1)+
  ylim(0.5, 1.1)+
  theme_publication+
  labs(title="EVmutation-+WT_CSTopCSTop", x="DTp56libD18-predictedFitness", y="DTp58libD21-predictedFitness")
ggsave(plot=a1, filename = "Sa888889909_EVmutationCTCT.svg", width=3, height=3)

#2) EVmutation, CLADETopTopTop
t<-SaCas9_888889909_prediction %>%   mutate(samplename=paste0(sg, "_", predictor, "_", sampling_r2) ) %>% 
  select(AACombo, sg, `InTrainingData?`, PredictedFitness, Fitness_n, rank_n, samplename, predictor, sampling_r1, sampling_r2) 
left_join(left_join(t %>% group_by(samplename) %>% mutate(prank_p=percent_rank(PredictedFitness), prank_n=percent_rank(Fitness_n)) %>% filter(prank_n>=0.95 & `InTrainingData?`=="YES")
                    %>% group_by(samplename) %>% summarise(n=n()), 
                    t %>% filter(`InTrainingData?`=="YES") %>% group_by(samplename) %>% summarise(s=n())) %>% mutate(n/s),
          t %>% group_by(samplename) %>% mutate(rank_p=rank(-PredictedFitness, ties.method  ="min")) %>%
            filter(rank_n<=50, rank_p<=50) %>%
            group_by(samplename) %>% summarise(t=n()))
t %>% filter(rank_n<=50, `InTrainingData?`=="YES") %>%
  group_by(samplename) %>% summarise(v=n())
### only ploting predicted fitness here
t<-t %>% filter(predictor=="EVmutation", sampling_r2=="CLADETopTopTop") %>%
  select(AACombo, PredictedFitness, predictor, sg,  sampling_r1, sampling_r2) %>%
  spread(sg, PredictedFitness)
t<-t%>% mutate(rank1=rank(-DTp56libD18, ties.method  ="min"),
               rank2=rank(-DTp58libD21, ties.method  ="min")) 
EVmutation_Top_CTTT<- t %>%   filter(rank1<=50 & rank2<=50) %>% select(AACombo) %>% unlist()
a1<-t %>% ggplot(aes(x=DTp56libD18, y=DTp58libD21, label=AACombo))+
  geom_point(alpha=0.05)+
  geom_point(data=t %>% filter(rank1<=50), colour="cyan", alpha=0.5)+
  geom_point(data=t %>% filter(rank2<=50), colour="deeppink", alpha=0.4)+
  geom_point(data=t %>% filter(AACombo=="NAL"), color="green")+
  geom_point(data= t %>% filter(AACombo %in% v_train$AACombo[v_train$sampling_r2=="CLADETopTopTop" & v_train$sg=="DTp56libD18" & v_train$predictor=="EVmutation"]), pch=1, colour="blue", size=2.5)+
  geom_point(data= t %>% filter(AACombo %in% v_train$AACombo[v_train$sampling_r2=="CLADETopTopTop" & v_train$sg=="DTp58libD21" & v_train$predictor=="EVmutation"]), pch=1, colour="red", size=1.5)+
  geom_text_repel(data= t %>% 
                    filter(rank1<=50 & rank2<=50), 
                  aes(label=AACombo), max.overlaps = 100, size=3, min.segment.length = 0.01)+
  geom_text_repel(data=t %>% filter(AACombo=="NAL"), color="green", size=3.5)+
  xlim(0.5, 1)+
  ylim(0.5, 1.1)+
  theme_publication+
  labs(title="EVmutation-+WT_CSTopTopTop", x="DTp56libD18-predictedFitness", y="DTp58libD21-predictedFitness")
ggsave(plot=a1, filename = "Sa888889909_EVmutationCTTT.svg", width=3, height=3)

### EGFP validation
egfp<-read_csv("SaCas9_888889909_validation_egfp.csv")
a1<-egfp %>% ggplot(aes(x=sg, y=EGFP,group=AACombo, fill=AACombo))+
  geom_bar(stat = "summary", fun.y = "mean", position = position_dodge(width = 0.8))+
  geom_point(position = position_dodge(width = 0.8))+
  scale_fill_manual(values = c("NAL"="gray70", "GAL"="skyblue"))+
  theme_publication+
  labs(x="sgRNA", y="EGFP disruption(%)")
ggsave(plot=a1, filename = "egfp.svg", width=3, height=3)