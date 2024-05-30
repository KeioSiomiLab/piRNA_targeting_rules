.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)


# haploid analysis
haploid_analysis <- function(keyname){
  repmask <- read_tsv(paste0("/media/pericles/TEfind/flam/hap/rmblast2/flam_",keyname,".fasta.out.gff"),
                      col_names = c("chr","source", "dis1","start","end", "score","strand","dis2","name"),
                      skip = 3)
  repmask2 <- repmask %>% 
    mutate(name2 = name %>% str_remove('Target "Motif:') %>% str_remove('"') %>% str_replace_all(" ","_")) %>% 
    select(chr,start,end,name2,score,strand) 
  repmask2 %>% 
    write_tsv(paste0("/media/pericles/TEfind/flam/hap/rmblast2/",keyname,"_TE.bed"), col_names = FALSE)
  repmask2 %>% filter(strand == "+") %>% 
    write_tsv(paste0("/media/pericles/TEfind/flam/hap/rmblast2/",keyname,"_TE_plus.bed"), col_names = FALSE)
  repmask2 %>%  filter(strand == "-") %>% 
    write_tsv(paste0("/media/pericles/TEfind/flam/hap/rmblast2/",keyname,"_TE_minus.bed"), col_names = FALSE)
  
  genes <- read_tsv(paste0("/media/pericles/TEfind/flam/hap/gene_res_",keyname,".txt"),
                    col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                  "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
  genes %>% 
    select(sseqid, sstart, send, qseqid) %>% 
    write_tsv(paste0("/media/pericles/TEfind/flam/hap/gene_res",keyname,".bed"),
              col_names = FALSE)
}
haploid_analysis(keyname = "hap1")
haploid_analysis(keyname = "hap2")
