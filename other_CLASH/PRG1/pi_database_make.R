#######################################################################################################
# piRNA database making
#######################################################################################################
#piRNA for c elegans
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/hermione/PRG1_CLASH/")
library(tidyverse)
piRNA1 <- read_tsv("Process/SRR2140770_comp.tsv",
                   col_names = c("piRNA", "seq")) %>% 
  mutate(sample = "piRNA1",
         length = str_length(seq))
piRNA2 <- read_tsv("Process/SRR538357_comp.tsv",
                   col_names = c("piRNA", "seq")) %>% 
  mutate(sample = "piRNA2",
         length = str_length(seq))

piRNA <- bind_rows(piRNA1,piRNA2) %>% 
  separate(piRNA, into = c("name","count"),sep ="_") %>% 
  mutate(count = count %>% as.integer()) %>% 
  group_by(seq) %>% summarise(count = sum(count)) %>% ungroup() %>% 
  filter(seq != "AAAAAAAAAAAAAAAAAAAAA")
piRNAtotal <- piRNA %>% arrange(desc(count)) %>% 
  mutate(piRNA = paste0("piRNA",row_number(),"_",count)) %>% 
  select(piRNA,seq) 
writepiRNA <- piRNAtotal %>% 
  mutate(annotation2 = paste0(">", piRNA), index = row_number()) %>% 
  select(annotation2,seq,index) %>% 
  pivot_longer(cols = -index, names_to = "type", values_to = "fasta") %>% 
  arrange(index,type) %>% select(fasta)
piRNAtotal %>% write_tsv("database/piRNA_all.tsv",col_names = FALSE)
writepiRNA %>% write_tsv("database/piRNA_all.fasta",col_names = FALSE)

#####################################################
#database for CLASH
filedir <- "/media/hermione/PRG1_CLASH/Process2/"

files3 <- dir(filedir, pattern = "_comp2.tsv$", full.names = TRUE)
files4 <- files3 %>%
  str_replace_all("_comp2.tsv$", "") %>% 
  str_replace(paste0(filedir, "/"), "") 
for (i in seq_along(files3)) {
  read_tsv(files3[[i]], col_names = c("name", "clashseq")) %>% 
    mutate(seq = clashseq %>% str_sub(1,21)) %>% left_join(piRNA,by="seq") %>% 
    drop_na()%>% mutate(seq2= clashseq %>% str_sub(start = 22),
                        left_length = seq2 %>% str_length()) %>% 
    rename(piRNAseq = seq) %>% 
    filter(left_length > 21) %>% 
    write_tsv(paste0("Process2/", files4[[i]], "_pimap.tsv")) %>% 
    mutate(annotation2 = paste0(">", name), index = row_number()) %>% 
    select(annotation2,seq2,index) %>% 
    pivot_longer(cols = -index, names_to = "type", values_to = "fasta") %>% 
    arrange(index,type) %>% select(fasta) %>%
    write_tsv(paste0("Process2/", files4[[i]], "_comp_tf.fasta"),col_names = FALSE)
}
