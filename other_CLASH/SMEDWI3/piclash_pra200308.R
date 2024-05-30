#######################################################################################################
# trim piRNA and generate target sequences
#######################################################################################################
# /usr/bin/Rscript  piclash_pra200308.R Process2/{/.}_map.sam Process2/{/.}_comp.tsv database/piRNA_all.tsv Process2/{/.}_comp_tf.fasta
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)

# hybdir <- "/media/hermione/SMEDWI3/Process2/SRR8842976_map.sam"
# CLASHdir <- "/media/hermione/SMEDWI3/Process2/SRR8842976_comp.tsv"
# piRNAdir <- "/media/hermione/SMEDWI3/database/piRNA_all.tsv"
# outdir <- "/media/hermione/SMEDWI3/Process2/SRR8842976_comp_tf.fasta"

hybdir <- commandArgs(trailingOnly = TRUE)[1]
CLASHdir <- commandArgs(trailingOnly = TRUE)[2]
piRNAdir <- commandArgs(trailingOnly = TRUE)[3]
outdir <- commandArgs(trailingOnly = TRUE)[4]

piRNA <- read_tsv(piRNAdir, col_names = c("name", "piseq"))

hyb1 <- read_tsv(hybdir, 
                 col_names = c("name", "flag", "target", "start", "qual", "cigar")) %>% 
  drop_na() %>% 
  filter(flag != 4 & flag !=16) %>% 
  left_join(piRNA, by="name")

hyb2 <- hyb1 %>% filter(start ==1) %>% 
  filter(!(str_detect(cigar, "I|D|S")))

CLASH <- read_tsv(CLASHdir, col_names = c("target", "seq"))
CLASH2 <- CLASH %>% mutate(totallen = seq %>% str_length())

Sub_CLASH <- hyb2 %>% left_join(CLASH, by = "target") %>% 
  mutate(match = cigar %>% str_replace("M", "") %>% as.integer(),
         new_seq = seq %>% str_sub(start=match+1),
         matchseq = seq %>% str_sub(end = match),
         new_length = new_seq %>% str_length(),
         direction = "pi_tag") %>% 
  filter(matchseq == piseq) %>% 
  filter(new_length >= 15) %>% 
  unite(name, target, direction, col="new_name",sep ="@") %>% 
  select(new_name,new_seq)

hyb3 <- hyb1 %>% left_join(CLASH2, by = "target") %>% 
  mutate(pilen = str_length(piseq)) %>% 
  filter(pilen + start -1 == totallen) %>% 
  filter(!(str_detect(cigar, "I|D|S")))

Sub_CLASH2 <- hyb3 %>% 
  mutate(match = cigar %>% str_replace("M", "") %>% as.integer(),
         new_seq = seq %>% str_sub(end = start - 1),
         matchseq = seq %>% str_sub(start = start),
         new_length = new_seq %>% str_length(),
         direction = "tag_pi") %>% 
  filter(matchseq == piseq) %>% 
  filter(new_length >= 15) %>% 
  unite(name, target, direction, col="new_name",sep ="@") %>% 
  select(new_name,new_seq)


fasta <- Sub_CLASH %>%  bind_rows(Sub_CLASH2) %>%
  mutate(annotation2 = paste0(">", new_name), index = row_number()) %>% 
  pivot_longer(new_seq:annotation2,  names_to = "type", values_to = "fasta") %>% 
  arrange(index,type) %>% select(fasta)
fasta %>% write_tsv(outdir, col_names = FALSE) 
