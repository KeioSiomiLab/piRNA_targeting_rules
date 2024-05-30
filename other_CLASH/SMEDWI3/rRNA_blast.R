
#######################################################################################################
# make rRNA site in genome
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)


test <- read_tsv("/media/hermione/SMEDWI3/database/raw/rRNA_align_out.txt",
                 col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                               "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>% 
  arrange(sseqid,sstart) %>% 
  mutate(strand = if_else(sstart > send, "antisense", "sense"),
         sstart2 = if_else(strand == "antisense", send, sstart),
         send2 = if_else(strand == "antisense", sstart, send)) %>% 
  mutate(strand2 = if_else(strand=="sense","+","-")) %>% 
  select(sseqid, sstart2, send2, strand2)%>% distinct() %>% 
  mutate(name = paste0("rRNA",row_number()),score = 0) %>% 
  select(sseqid, sstart2, send2, name,score,strand2) 

test %>% 
  write_tsv(paste0("/media/hermione/SMEDWI3/", "database/raw/rRNA_align.bed"),col_names = FALSE)


