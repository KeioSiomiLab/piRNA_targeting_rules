##########################################
# cage annotation
##########################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

library(tidyverse)

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

figoutdir <- "/media/hermione/CAGE/figures/" 
ctssdir <- "/media/hermione/CAGE/ctss/"
cagedir <- "/media/hermione/CAGE/genome/"
cage_anno <- read_tsv("/media/hermione/CAGE/figures/CAGE_annotate.bed",
         col_names = c("chr","start","end","remove1","score", "direction",
                       "remove2","remove3","remove4","anno","remove","remove5","distance")) %>%
  select(-contains("remove")) %>% mutate(type = "close")

cage_inter <- read_tsv("/media/hermione/CAGE/figures/CAGE_inter.bed",
                       col_names = c("chr","start","end","remove1","score", "direction",
                                     "remove2","remove3","remove4","anno","remove","remove5")) %>% 
  select(-contains("remove")) %>% mutate(type = "5UTR")

cage_remove <- cage_anno %>% filter(distance > 100) %>% 
  anti_join(cage_inter, by = c("chr","start","end","direction","anno"))

cage_anno2 <- bind_rows(cage_anno, cage_inter) %>%
  anti_join(cage_remove, by = c("chr","start","end","direction","anno")) %>% 
  select(-distance, -type) %>% distinct()

tss_anno <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_tss.tsv") %>% rename(genename = name) %>% 
  select(genename,transcript, genesymbol, tfsymbol)


cage_anno3 <- cage_anno2 %>%
  left_join(tss_anno %>% rename(anno = transcript), by = "anno") %>%
  group_by(chr,start,end, direction,genename,genesymbol) %>%
  summarise(anno = paste0(anno,collapse = "&"), tfsymbol = paste0(tfsymbol,collapse = "&")) %>%
  ungroup()

cage_anno3 %>%
  write_tsv(paste0(figoutdir,"ctss_annotation.tsv"))
  
# cage data for SQANTI3
all_peaks <- read_tsv("/media/hermione/CAGE/figures/all_peaks_merged.bed",
                      col_names = FALSE)

all_peaks %>% mutate(name = paste0("mergedpeak",row_number()), strand = X4 %>% str_sub(start = 1L, end = 1L),score = 0) %>% 
  select(X1, X2, X3, name,score, strand) %>%  
  write_tsv("/media/hermione/CAGE/figures/All.samples.tagClusters_SQANTI.bed", col_names = FALSE,escape = "none")
  
EGFP_peaks <- read_tsv("/media/hermione/CAGE/figures/siEGFP.tagClusters.qLow0.1_qUp0.9.bed",
                      skip = 1,col_names = FALSE)
EGFP_peaks %>% mutate(name = paste0(X1, ":",X2,"..",X3,",",X6)) %>% 
  select(X1, X2, X3, name,X5, X6,X7,X8) %>% 
  write_tsv("/media/hermione/CAGE/figures/siEGFP.tagClusters_SQANTI.bed", col_names = FALSE,escape = "none")

Piwi_peaks <- read_tsv("/media/hermione/CAGE/figures/siPiwi.tagClusters.qLow0.1_qUp0.9.bed",
                       skip = 1,col_names = FALSE)
Piwi_peaks %>% mutate(name = paste0(X1, ":",X2,"..",X3,",",X6)) %>% 
  select(X1, X2, X3, name,X5, X6,X7,X8) %>% 
  write_tsv("/media/hermione/CAGE/figures/siPiwi.tagClusters_SQANTI.bed", col_names = FALSE,escape = "none")
