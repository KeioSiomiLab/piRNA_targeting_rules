

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))

outdir <- "/media/hermione/piRNAmodel/"
outname <- "/media/hermione/piRNAmodel/vienna/"
#######################################################################################################
# processing piRNA mapping to the models
#######################################################################################################
piRNA_seq <- read_tsv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/piRNA/", "piRNA_CPM.tsv")) %>% 
  select(name,seq, pilen,CPM_avg,type3)

TE_rnaplex3 <- read_tsv(paste0(outname,"model_RNAplex","_vienna.tsv"))
gene_rnaplex3 <- read_tsv(paste0(outname,"model_RNAplex2","_vienna.tsv"))

gene_remove <- read_tsv(paste0(outdir ,"rnaplex/","gene_remove1_overlap.tsv"),
                        col_names = c("dis1","dis2","dis3","readname","targetname","dis4","dis5","TE","annotype")) %>% 
  filter(str_detect(annotype,"TE")) %>% select(-contains("dis"))

All_rnaplex3 <- bind_rows(gene_rnaplex3 %>% anti_join(gene_remove, by = "readname") %>% mutate(pairtype = "gene"),
                          TE_rnaplex3 %>% mutate(pairtype = "TE")) %>% 
  left_join(piRNA_seq %>% rename(piRNAname = name) %>% select(piRNAname, type3), by = "piRNAname")


All_rnaplex3 %>% 
  write_tsv(paste0(outname,"All_model_vienna.tsv"))

