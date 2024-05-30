#######################################################################################################
# filtering pair and preparing for RNAplex
#######################################################################################################
# /usr/bin/Rscript  miclash_tf_sep3.R Process3/{/.}_tf_map.sam Process3/{/.}_comp_tf.tsv ref/hOH7-microRNA.tsv Process3/{/.}_chimera
# hybdir <- "/media/hermione/AGO1_CLASH/Process3/SRR959751_tf_map.sam"
# CLASHdir <- "/media/hermione/AGO1_CLASH/Process3/SRR959751_comp_tf.tsv"
# miRNAdir <- "/media/hermione/AGO1_CLASH/ref/hOH7-microRNA.tsv"
# outname <- "/media/hermione/AGO1_CLASH/Process3/SRR959751_chimera"


.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)

hybdir <- commandArgs(trailingOnly = TRUE)[1]
CLASHdir <- commandArgs(trailingOnly = TRUE)[2]
miRNAdir <- commandArgs(trailingOnly = TRUE)[3]
outname <- commandArgs(trailingOnly = TRUE)[4]

hyb1 <- read_tsv(hybdir, 
                 col_names = c("target", "flag", "target2", "start", "qual", "cigar")) %>% 
  drop_na() %>% 
  filter(flag != 4 & flag !=16) 

miRNA <- read_tsv(miRNAdir, col_names = c("miRNA", "miRNAseq"))

CLASH <- read_tsv(CLASHdir, col_names = c("target", "targetseq"))

hyb_all_pre <- hyb1 %>% left_join(CLASH, by="target") %>% 
  separate(target, sep = "@",into = c("miRNA", "readname", "direction")) %>% 
  left_join(miRNA, by = "miRNA")
cigar_index <- hyb_all_pre %>% select(cigar) %>% distinct() %>% 
  mutate(pack = map(cigar,function(x){str_split(x,pattern = "[:upper:]")}),
         titania = map(cigar,function(x){str_split(x,pattern = "[:digit:]")})) %>%
  unnest_legacy() %>% 
  mutate(pack2 = map(pack, function(x) {x[which(x != "")]}),
         titania2 = map(titania, function(x) {x[which(x != "")]}),
         titania3 = map(titania2, function(x) {(which(x=="I"|x=="S"))}),
         end = map2(pack2, titania3, function(x,y) {replace(x,y,0) %>% as.integer() %>% sum()})) %>% 
  select(-contains("pack"),-contains("titania")) %>% 
  unnest_legacy()
hyb_all <- hyb_all_pre %>% left_join(cigar_index,by = "cigar") %>% 
  mutate(milen = str_length(miRNAseq),
         bedstart = start-1,
         bedend =start + end-1,
         target2 = str_remove_all(target2, "\\-"),
         targetlen = str_length(targetseq)) %>% 
#  filter(targetlen > milen) %>% 
  select(readname,miRNA,milen,miRNAseq,target2,bedstart,bedend,targetlen,targetseq,cigar,direction)


hyb_all %>% write_tsv(paste0(outname, ".tsv"))

hyball2 <- hyb_all %>% select(readname,miRNA,miRNAseq,targetseq) %>% 
  distinct()

hyball3 <- hyball2 %>% mutate(index = row_number(),miRNA = paste0(">",miRNA), readname = paste0(">",readname)) %>% 
  pivot_longer(cols = c("miRNA","miRNAseq","readname","targetseq"), names_to = "type", values_to = "fasta") %>% 
  select(fasta)
hyball3 %>% write_tsv(paste0(outname, "_rnaplex.fasta"), col_names = FALSE)


