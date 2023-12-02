library(ggplot2)
library(tidyverse)
library(ggpmisc)
library(ggseqlogo)
library(readxl)
library(Cairo)
library(stringdist)
library(ggsci)
theme_publication<-theme_classic()+
  theme(text = element_text(size=12, family = "Arial"),
        panel.background = element_rect(colour="black", size=0.5))

Cas9_xl<-read_excel("41592_2019_473_MOESM4_ESM_CombiSEAL_supp_Fig2.xlsx")
WT="RQKETQKR"
non_target_order=c("---", "-R-", "-H-", "--A", "-RA", "-HA", "-AA", "A--", "AR-", "AH-", "A-A", "ARA", "AHA", "AAA")
target_order=c( "-----", "A----", "-A---", "AA---", "----A", "A---A", "-A--A", "AA--A", "--QGP", 
                "A-QGP", "-AQGP", "AAQGP", "--VRE", "A-VRE", "-AVRE", "AAVRE", "--AWE", "A-AWE", "-AAWE" ,"AAAWE",
                "--MVA", "A-MVA", "-AMVA", "AAMVA", "--KSA", "A-KSA", "-AKSA", "AAKSA", "--RK-", "A-RK-", "-ARK-",  "AARK-", 
                "--CRE", "A-CRE" , "-ACRE", "AACRE", "--QW-", "A-QW-", "-AQW-", "AAQW-", "--LGA", "A-LGA", "-ALGA", "AALGA",
                "--WDE", "A-WDE" , "-AWDE" , "AAWDE", "--HL-", "A-HL-", "-AHL-", "AAHL-", "--VWA", "A-VWA", "-AVWA", "AAVWA",
                "--RRA" , "A-RRA", "-ARRA", "AARRA", "--GDE", "A-GDE", "-AGDE", "AAGDE" , "--MRA", "A-MRA",  "-AMRA", "AAMRA")
colnames(Cas9_xl)<-Cas9_xl[4,]
Cas9_xl<-Cas9_xl[5:956, ]
Cas9_xl[, c("RFPsg5 ON", "RFPsg5 OFF5-2", "RFPsg8 ON", "RFPsg8 OFF5")]<-sapply(Cas9_xl[, c("RFPsg5 ON", "RFPsg5 OFF5-2", "RFPsg8 ON", "RFPsg8 OFF5")], as.numeric)
sg5_Fit<-Cas9_xl %>% mutate(AACombo=paste0(`661`, `695`, `848`, `923`, `924`, `926`, `1003`, `1060`, sep="")) %>% select(AACombo, `RFPsg5 ON`)
colnames(sg5_Fit)<-c("AACombo", "Fitness")
sg8_Fit<-Cas9_xl %>% mutate(AACombo=paste0(`661`, `695`, `848`, `923`, `924`, `926`, `1003`, `1060`, sep="")) %>% select(AACombo, `RFPsg8 ON`)
colnames(sg8_Fit)<-c("AACombo", "Fitness")

### make mutation list for EvoEF
mut_nontation= Cas9_xl %>% mutate(mut_no=paste0("RB661", `661`, ",", "QB695", `695`, ",", "KB848", `848`, ",",
                           "EB923", `923`, ",", "TB924", `924`, ",", "QB926", `926`, ",",
                            "KB1003", `1003`, ",", "RB1060", `1060`, ";")) %>%select(mut_no)
write_csv(mut_nontation, "/media/achu/新增磁碟區/MLDE/2021_11_09_new_MLDE/SpCas9_586-1255/SpCas9_mutlist.txt", col_names = F, quote = "none")


log_text<-read_delim("/media/achu/新增磁碟區/MLDE/2021_11_09_new_MLDE/SpCas9_586-1255/SpCas9_log.txt", col_names = F)
log_text<-log_text[, 1:2]
colnames(log_text)<-c("PDB", "DG")
WT_DG<-  log_text[grepl("WT", log_text$PDB), ] %>% summarise(mean(DG)) %>% unlist()
log_text<-log_text[!grepl("WT", log_text$PDB), ]
Cas9_xl <- Cas9_xl %>% mutate(AACombo=paste0(`661`, `695`, `848`, `923`, `924`, `926`, `1003`, `1060`, sep=""))
log_text<-bind_cols(Cas9_xl[, "AACombo"], log_text)
log_text<-log_text %>% mutate(DDG_stab=DG - WT_DG)
write_csv(log_text, "SpCas9_R661Q695K848E923T924Q926K1003R1060_EvoEF.csv")

fit_tab<-Cas9_xl[, c("AACombo", "RFPsg5 ON", "RFPsg8 ON")]
colnames(fit_tab)<-c("AACombo", "Sg5", "Sg8")
### Do min Mac
fit_tab$Sg8[fit_tab$Sg8< -2 & !is.na(fit_tab$Sg8)]<- -2
fit_tab$Sg5<- (fit_tab$Sg5 - min(fit_tab$Sg5, na.rm = T))/(max(fit_tab$Sg5, na.rm = T) - min(fit_tab$Sg5, na.rm = T))
fit_tab$Sg8<- (fit_tab$Sg8 - min(fit_tab$Sg8, na.rm = T))/(max(fit_tab$Sg8, na.rm = T) - min(fit_tab$Sg8, na.rm = T))

fit_tab %>% ggplot(aes(x=Sg5))+geom_histogram()
fit_tab %>% ggplot(aes(x=Sg8))+geom_histogram()
write_csv(fit_tab, "SpCas9_fitness.csv")
