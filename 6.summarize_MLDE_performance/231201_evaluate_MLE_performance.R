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

#### SaCas9_WED
setwd("data/MLDE_train_test_datasets")
Sa_Fit<-read_csv("SaCas9_8888889909_Fitness_norm.csv")
test_sample<-read_csv("Sa8000_test_sample2.data", col_names = F)
Sa_test_sample<-unique(unlist(test_sample[, 1]))

#### SaCas9_WED-PI
Sa_Fit_DT<-read_csv("SaCas9_887888887985986988989991_Fitness_norm.csv")
test_sample_DT<-read_csv("SaCas9_887888889985986988989991_fixed_test_set.txt", col_names = F)
test_sample_DT<-unlist(test_sample_DT[, 1])

#### SpCas9 
Sp_Fit<-read_csv("SpCas9_fitness_norm.csv")
### import test samples
Sp_test_sample<-read_csv("SpCas9_R661Q695K848E923T924Q926K1003R1060_fixed_test_set.txt", col_names = F)
Sp_test_sample<-unlist(Sp_test_sample[, 1])

### PhoQ ###
Fit<-read_csv("PhoQ_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
PhoQ_Fit<-Fit
test_sample<-read_csv("PhoQ_test_sample.data")
colnames(test_sample)<-c("AACombo", "Fitness", "Fitness_n")
PhoQ_test_sample<-test_sample %>% select(AACombo) %>% unlist()

### GFP ###
Fit<-read_csv("GFP_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
GFP_Fit<-Fit
GFP_test_sample<-read_csv("GFP_test_sample.data") %>% select(AACombo) %>% unlist()

### PABP ###
Fit<-read_csv("PABP_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
PABP_Fit<-Fit
PABP_test_sample<-read_csv("PABP_test_sample.data") %>% select(AACombo) %>% unlist()

### UBE ###
Fit<-read_csv("UBE_zeroshot.csv")
Fit<-Fit[, 2:8]
colnames(Fit)[1]<-"AACombo"
UBE_Fit<-Fit
UBE_test_sample<-read_csv("UBE_test_sample.data")%>% select(AACombo) %>% unlist()


#####
compute_performance<-function(samplenames, Fit, test_sample){
  training_variants<-c()
  max_pre<-c()
  cors<-c()
  t5<-c()
  NDCG<-c()
  v_t5<-c()
  total_variants<-c()
  top_50<-c()
  top_20<-c()
  top_1_percent<-c()
  top_1_percent_vno<-c()
  top_50_training<-c()
  top_20_training<-c()

  for (i in samplenames){
    path<-list.dirs(path = i, recursive = F)
    print(path)
    if (length(path)==0){
    p<-read_csv(paste0(i, "/PredictedFitness.csv"))
    } else if (length(path)==1){
    p<-read_csv(paste0(path, "/PredictedFitness.csv"))}
    else {
      p<-read_csv(paste0(path[length(path)], "/PredictedFitness.csv"))}
      p<-left_join(Fit, p)
      
      ### remove missing data rows 
      p<- p %>% filter(!is.na(Fitness_n))
      total_variants<-c(total_variants, p %>% nrow())
      training_variants<-c(training_variants, p %>% filter(`InTrainingData?`=="YES") %>% nrow())
      v_t5<-c(v_t5, p%>% mutate(rank_n=percent_rank(Fitness_n)) %>%
                filter(rank_n>=0.95, `InTrainingData?`=="YES") %>% nrow() / 
                p%>% mutate(rank_n=percent_rank(Fitness_n)) %>%
                filter(`InTrainingData?`=="YES") %>% nrow())
      top_1_percent<-c(top_1_percent, p %>% mutate(rank_n=percent_rank(Fitness_n), rank_p=percent_rank(PredictedFitness)) %>% filter(rank_n>=0.99 & rank_p>=0.99) %>% nrow())
      top_1_percent_vno<-c(top_1_percent_vno, p %>% mutate(rank_n=percent_rank(Fitness_n)) %>% filter(rank_n>=0.99) %>% nrow())
      top_50<-c(top_50, p %>% mutate(rank_n=rank(-Fitness_n, na.last = T, ties.method = "min"),
                                     rank_p=rank(-PredictedFitness, na.last = T,ties.method = "min"))
                %>% filter(rank_n<=50, rank_p<=50) %>% nrow())
      top_20<-c(top_20, p %>% mutate(rank_n=rank(-Fitness_n, na.last = T, ties.method = "min"),
                                     rank_p=rank(-PredictedFitness, na.last = T,ties.method = "min"))
                %>% filter(rank_n<=20, rank_p<=20) %>% nrow())
      top_50_training<-c(top_50_training, p %>% mutate(rank_n=rank(-Fitness_n, na.last = T, ties.method = "min"),
                                                       rank_p=rank(-PredictedFitness, na.last = T,ties.method = "min"))
                         %>% filter(rank_n<=50, rank_p<=50) %>%filter(`InTrainingData?`=="YES") %>% nrow())
      top_20_training<-c(top_20_training,  p %>% mutate(rank_n=rank(-Fitness_n, na.last = T, ties.method = "min"),
                                                        rank_p=rank(-PredictedFitness, na.last = T,ties.method = "min"))
                         %>% filter(rank_n<=20, rank_p<=20) %>%filter(`InTrainingData?`=="YES") %>% nrow())
      ### remove points present in training data ###
    p<- p %>% filter(AACombo %in% test_sample)
    
    p<- p %>% filter(`InTrainingData?`=="NO")
    max_pre<-c(max_pre, max(p$PredictedFitness))
    ### calculate correlation
    r<-cor.test(p$PredictedFitness, p$Fitness_n, method = "spearman")
    cors<-c(cors, r$estimate)
    ### calculate NDCG
    ndcg<-p %>% mutate(rank_p=rank(-PredictedFitness), rank_n=rank(-Fitness_n)) %>%
      mutate(pDGC= Fitness_n/log(rank_p+1), IDGC = Fitness_n/log(rank_n +1)) %>%
      summarise(pDGC=sum(pDGC, na.rm = T), IDGC=sum(IDGC, na.rm = T)) %>%
      mutate(NDGC=pDGC/IDGC) %>% select(NDGC) %>% ungroup() %>% unlist()
    NDCG<-c(NDCG, ndcg)
    ### calculate enrichment
    o<-p %>% filter(!is.na(Fitness_n))%>% mutate(rank_n=percent_rank(Fitness_n), rank_p=percent_rank(PredictedFitness)) %>%
      filter(rank_n>=0.95 & rank_p>=0.95) %>% nrow()
    n<-p %>% filter(!is.na(Fitness_n)) %>% nrow
    t5<-c(t5, o*400/n)
  }
  stats_tab<-tibble(samplenames=samplenames, cors, t5, NDCG, max_pre, v_t5, top_1_percent, top_50, top_20, top_50_training, top_20_training, training_variants, total_variants, top_1_percent_vno)
  return(stats_tab)}

setwd("data/execute_MLDE_multi_round/direct_ranking_fitness_execute_mlde")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))
write_csv(stats_tab, paste0(Sys.Date(), "_zeroshot_stats_tab.csv"))



setwd("data/execute_MLDE_multi_round/execute_mlde_zeroshot")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))
write_csv(stats_tab, paste0(Sys.Date(), "_zeroshot_stats_tab.csv"))


setwd("data/execute_MLDE_multi_round/execute_mlde_zeroshot_WT")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))
write_csv(stats_tab, paste0(Sys.Date(), "_zeroshot_stats_tab.csv"))

setwd("/media/achu/新增磁碟區/CLADE/2nd_round_execute_mlde_WT")

write_csv(stats_tab, paste0(Sys.Date(), "_2nd_round_zeroshot_stats_tab.csv"))

setwd("data/execute_MLDE_multi_round/2nd_round_execute_mlde")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))
write_csv(stats_tab, paste0(Sys.Date(), "_2nd_round_zeroshot_stats_tab.csv"))

setwd("data/execute_MLDE_multi_round/3rd_round_execute_mlde")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))
write_csv(stats_tab, paste0(Sys.Date(), "_3rd_round_zeroshot_stats_tab.csv"))

setwd("data/execute_MLDE_multi_round/3rd_round_execute_mlde_WT")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))
write_csv(stats_tab, paste0(Sys.Date(), "_3rd_round_zeroshot_stats_tab.csv"))

setwd("data/execute_MLDE_multi_round/4th_round_execute_mlde_WT")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))
write_csv(stats_tab, paste0(Sys.Date(), "_4th_round_zeroshot_stats_tab.csv"))

setwd("data/execute_MLDE_multi_round/4th_round_execute_mlde")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))
write_csv(stats_tab, paste0(Sys.Date(), "_4th_round_zeroshot_stats_tab.csv"))



setwd("data/execute_MLDE_benchmark/execute_1st_round_MLDE_benchmark")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("SpCas9_Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("SpCas9_Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("DTp56libD18", samplenames)]
Fit_temp<-Sa_Fit %>% select(AACombo, Fitness_n=DTp56libD18)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sa_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("DTp58libD21", samplenames)]
Fit_temp<-Sa_Fit %>% select(AACombo, Fitness_n=DTp58libD21)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sa_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("PhoQ", samplenames)]
Fit_temp<-PhoQ_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, PhoQ_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GFP", samplenames)]
Fit_temp<-GFP_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, GFP_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("UBE", samplenames)]
Fit_temp<-UBE_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, UBE_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("PABP", samplenames)]
Fit_temp<-PABP_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, PABP_test_sample))

write_csv(stats_tab, paste0(Sys.Date(), "_zeroshot_stats_tab.csv"))

setwd("data/execute_MLDE_benchmark/execute_4th_round_MLDE_benchmark")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("SpCas9_Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("SpCas9_Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("DTp56libD18", samplenames)]
Fit_temp<-Sa_Fit %>% select(AACombo, Fitness_n=DTp56libD18)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sa_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("DTp58libD21", samplenames)]
Fit_temp<-Sa_Fit %>% select(AACombo, Fitness_n=DTp58libD21)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sa_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("PhoQ", samplenames)]
Fit_temp<-PhoQ_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, PhoQ_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GFP", samplenames)]
Fit_temp<-GFP_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, GFP_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("UBE", samplenames)]
Fit_temp<-UBE_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, UBE_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("PABP", samplenames)]
Fit_temp<-PABP_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, PABP_test_sample))
write_csv(stats_tab, paste0(Sys.Date(), "_zeroshot_stats_tab.csv"))


setwd("data/execute_MLDE_benchmark/execute_2nd_round_MLDE_benchmark")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("SpCas9_Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("SpCas9_Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("DTp56libD18", samplenames)]
Fit_temp<-Sa_Fit %>% select(AACombo, Fitness_n=DTp56libD18)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sa_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("DTp58libD21", samplenames)]
Fit_temp<-Sa_Fit %>% select(AACombo, Fitness_n=DTp58libD21)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sa_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("PhoQ", samplenames)]
Fit_temp<-PhoQ_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, PhoQ_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GFP", samplenames)]
Fit_temp<-GFP_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, GFP_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("UBE", samplenames)]
Fit_temp<-UBE_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, UBE_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("PABP", samplenames)]
Fit_temp<-PABP_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, PABP_test_sample))

write_csv(stats_tab, paste0(Sys.Date(), "_zeroshot_stats_tab.csv"))

setwd("data/execute_MLDE_benchmark/execute_3rd_round_MLDE_benchmark")
samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("SpCas9_Sg5", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg5)
stats_tab<-compute_performance(samplenames, Fit_temp, Sp_test_sample)

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("SpCas9_Sg8", samplenames)]
Fit_temp<-Sp_Fit %>% select(AACombo, Fitness_n=Sg8)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sp_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("ACAAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=ACAAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GGTGGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=GGTGGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("TGGAGT", samplenames)]
Fit_temp<-Sa_Fit_DT %>% select(AACombo, Fitness_n=TGGAGT)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, test_sample_DT))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("DTp56libD18", samplenames)]
Fit_temp<-Sa_Fit %>% select(AACombo, Fitness_n=DTp56libD18)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sa_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("DTp58libD21", samplenames)]
Fit_temp<-Sa_Fit %>% select(AACombo, Fitness_n=DTp58libD21)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, Sa_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("PhoQ", samplenames)]
Fit_temp<-PhoQ_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, PhoQ_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("GFP", samplenames)]
Fit_temp<-GFP_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, GFP_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("UBE", samplenames)]
Fit_temp<-UBE_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, UBE_test_sample))

samplenames<-list.dirs(recursive = F)
samplenames<-samplenames[grepl("PABP", samplenames)]
Fit_temp<-PABP_Fit %>% select(AACombo, Fitness_n)
stats_tab<-bind_rows(stats_tab, compute_performance(samplenames, Fit_temp, PABP_test_sample))
write_csv(stats_tab, paste0(Sys.Date(), "_zeroshot_stats_tab.csv"))
