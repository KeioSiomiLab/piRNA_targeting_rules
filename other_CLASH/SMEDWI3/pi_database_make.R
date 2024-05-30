#######################################################################################################
# piRNA database
#######################################################################################################
# piRNA for planaria
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/hermione/SMEDWI3/")
library(tidyverse)
piRNA1 <- read_tsv("Process/SRR8162647_comp.tsv",
                   col_names = c("piRNA", "seq")) %>% 
  mutate(sample = "piRNA1", length = str_length(seq)) %>%
  filter(length <= 35) %>% filter(length>=26) %>% 
  mutate(count = piRNA %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer(),
         CPM = 1000000*count/sum(count)) %>% filter(count >= 3) %>% filter(CPM >= 0.1)
piRNA2 <- read_tsv("Process/SRR8162648_comp.tsv",
                   col_names = c("piRNA", "seq")) %>% 
  mutate(sample = "piRNA2", length = str_length(seq)) %>%
  filter(length <= 35) %>% filter(length>=26) %>% 
  mutate(count = piRNA %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer(),
         CPM = 1000000*count/sum(count)) %>% filter(count >= 3) %>% filter(CPM >= 0.1)
piRNA3 <- read_tsv("Process/SRR8162649_comp.tsv",
                   col_names = c("piRNA", "seq")) %>% 
  mutate(sample = "piRNA3", length = str_length(seq)) %>%
  filter(length <= 35) %>% filter(length>=26) %>% 
  mutate(count = piRNA %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer(),
         CPM = 1000000*count/sum(count)) %>% filter(count >= 3) %>% filter(CPM >= 0.1)



piRNA <- bind_rows(piRNA1,piRNA2,piRNA3) %>% 
  select(seq, sample,CPM) %>%
  group_by(seq) %>%
  summarize(samplecount = length(sample),
            type = list(sample),
            CPM = mean(CPM)) %>% ungroup() %>% 
  filter(samplecount >= 2) %>%
  filter(seq != "AAAAAAAAAAAAAAAAAAAAA")
piRNAtotal <- piRNA %>% arrange(desc(CPM)) %>% 
  mutate(piRNA = paste0("piRNA",row_number())) %>% 
  select(piRNA,seq) 
writepiRNA <- piRNAtotal %>% 
  mutate(annotation2 = paste0(">", piRNA), index = row_number()) %>% 
  select(annotation2,seq,index) %>% 
  pivot_longer(cols = -index, names_to = "type", values_to = "fasta") %>% 
  arrange(index,type) %>% select(fasta)
piRNAtotal %>% write_tsv("database/piRNA_all.tsv",col_names = FALSE)
writepiRNA %>% write_tsv("database/piRNA_all.fasta",col_names = FALSE)

