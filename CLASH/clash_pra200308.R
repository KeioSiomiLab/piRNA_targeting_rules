# clash_pra200308
# /usr/bin/Rscript --silent --slave --vanilla clash_pra200308.R Process2/{/.}_map.sam Process2/{/.}_comp.tsv piRNA/anno/piRNA_all2.tsv Process2/{/.}_comp_tf.fasta
#######################################################################################################
# trim piRNa and generate target sequences
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

#hybdir <- "/media/pericles/CLASH/Process2/rep1_map.sam"
#CLASHdir <- "/media/pericles/CLASH/Process2/rep1_comp.tsv"
#piRNAdir <- "/media/pericles/CLASH/piRNA/anno/piRNA_all2.tsv"
#outdir <- "/media/pericles/CLASH/Process2/rep1_comp_tf.fasta"
  
suppressMessages(suppressWarnings(require(tidyverse)))

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

fasta <- Sub_CLASH %>% bind_rows(Sub_CLASH2) %>%
  mutate(annotation2 = paste0(">", new_name), index = row_number()) %>% 
  pivot_longer(new_seq:annotation2,  names_to = "type", values_to = "fasta") %>% 
  arrange(index,type) %>% select(fasta)
fasta %>% write_tsv(outdir, col_names = FALSE) 
