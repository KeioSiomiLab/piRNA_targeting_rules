.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)
#flam


namelist <- c("flam_hap1", "flam_hap2", "20A_hap1", "20A_hap2")
for (i in seq_along(namelist)){
  repmask <- read_tsv(paste0("/media/pericles/TEfind/flam/rmblast/",namelist[[i]],"_region.fasta.out.gff"),
                      col_names = c("chr","source", "dis1","start","end", "score","strand","dis2","name"),
                      skip = 3)
  repmask2 <- repmask %>% 
    mutate(name2 = name %>% str_remove('Target "Motif:') %>% str_remove('"') %>% str_replace_all(" ","_")) %>% 
    select(chr,start,end,name2,score,strand) 
  repmask2 %>% 
    write_tsv(paste0("/media/pericles/TEfind/flam/rmblast/",namelist[[i]],"_TE.bed"), col_names = FALSE)
  repmask2 %>% filter(strand == "+") %>% 
    write_tsv(paste0("/media/pericles/TEfind/flam/rmblast/",namelist[[i]],"_TE_plus.bed"), col_names = FALSE)
  repmask2 %>%  filter(strand == "-") %>% 
    write_tsv(paste0("/media/pericles/TEfind/flam/rmblast/",namelist[[i]],"_TE_minus.bed"), col_names = FALSE)
}

