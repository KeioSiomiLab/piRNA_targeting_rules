.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)

rmblast <- read_tsv("/media/pericles/TEfind/rmblast/dm6_chr.fasta.out.gff", skip = 3,
                       col_names = c("chr","source", "dis1","start","end", "score","strand","dis2","name"))
rmblast2 <- rmblast %>% separate(name,into = c("dis3","name","TEstart","TEend"), sep = " ") %>% 
  mutate(name = name %>% str_remove_all('"|Motif:'))
rmblast2 %>% select(chr, start, end, name, score, strand) %>% 
  write_tsv("/media/pericles/TEfind/rmblast/dm6_rmblast_mask.bed",col_names = FALSE)
rmblast2 %>% select(chr, start, end, name, score, strand, TEstart, TEend) %>% 
  write_tsv("/media/pericles/TEfind/rmblast/dm6_rmblast_mask.tsv")

#############################################

