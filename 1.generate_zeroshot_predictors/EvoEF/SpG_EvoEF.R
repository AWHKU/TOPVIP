### try the SpG dataset D1135L/S1136W/G1218K/E1219Q/R1335Q/T1337R
library(tidyverse)
library(readxl)
SpG<-read_delim("SpG_Cas9_FigS4.csv")
### fitness:row2-65, PAM=NGAT, NGCC, NGGG, NGTA) 
a<-SpG[2:65,c(1, 2, 3)]
b<-SpG[2:65,c(1, 2, 4)]
c<-SpG[2:65,c(1, 2, 5)]
colnames(a)<-c("Name", "AACombo", "Fitness")
colnames(b)<-c("Name", "AACombo", "Fitness")
colnames(c)<-c("Name", "AACombo", "Fitness")
NGS_tab<-bind_rows(a %>% mutate(rep="rep1"),
                   b %>% mutate(rep="rep2"),
                   c %>% mutate(rep="rep3")) 

### look at sequences ###
p<-a %>% mutate(AA=AACombo) %>% 
  separate(AA, c("p", "1135", "1136", "1218", "1219", "1335", "1337"), sep = "")
AA_list<-apply(p[, c("1135","1136", "1218", "1219", "1335","1337")], 2 , unique)
#$`1135`
#[1] "D" "V" "A" "L" "F" "W"
#$`1136`
#[1] "S" "A" "V" "L" "F" "W"
#$`1218`
#[1] "G" "R" "K" "S"
#$`1219`
#[1] "E" "H" "Q" "S" "V"
#$`1335`
#[1] "R" "Q"
#$`1337`
#[1] "T" "R" "K"
library(seqinr)
Cas9<-read.fasta(file="spCas9_uniprot.fasta", 
                 seqtype = c("AA"), as.string = T, seqonly = T)

### try the SpG dataset D1135L/S1136W/G1218K/E1219Q/R1335Q/T1337R

### make EvoEF mutlist
evoEF_mutlist<-apply(p[, c("1135", "1136", "1218", "1219", "1335", "1337")], 
                     1,
                      function(x) paste0("DB1135", x[1], ",",
                                         "SB1136", x[2], ",",
                                         "GB1218", x[3], ",",
                                         "EB1219", x[4], ",",
                                         "RB1335", x[5], ",",
                                         "TB1337", x[6], ";"))
write_lines(evoEF_mutlist, "SpG_evoEF_mutlist.txt")

log_text<-read_delim("SpG_log.txt", col_names = F)
log_text<-log_text[, 1:2]
colnames(log_text)<-c("PDB", "DG")
WT_DG<-  log_text[grepl("WT", log_text$PDB), ] %>% summarise(mean(DG)) %>% unlist()
log_text<-log_text[!grepl("WT", log_text$PDB), ]
log_text<-bind_cols(p[, "AACombo"], log_text)
log_text<-log_text %>% mutate(DDG_stab=DG - WT_DG)
write_csv(log_text, "SpG_evoEF_zeroshot.txt")