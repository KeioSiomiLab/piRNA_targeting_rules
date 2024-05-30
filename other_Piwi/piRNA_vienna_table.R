

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))

#outdir <- "/media/hermione/otherdata_compare/Aub_testis/"

outdir <- commandArgs(trailingOnly = TRUE)[1]

#######################################################################################################
# processing piRNA mapping to the models
#######################################################################################################
piRNA_seq <- read_tsv(paste0(outdir, "piRNA_CPM.tsv")) %>% 
  select(name,seq, pilen,CPM_avg)

TE_rnaplex3 <- read_tsv(paste0(outdir,"vienna/","model_RNAplex","_vienna.tsv"))
gene_rnaplex3 <- read_tsv(paste0(outdir,"vienna/","model_RNAplex2","_vienna.tsv"))

gene_remove <- read_tsv(paste0(outdir ,"rnaplex/","gene_remove1_overlap.tsv"),
                        col_names = c("dis1","dis2","dis3","readname","targetname","dis4","dis5","TE","annotype")) %>% 
  filter(str_detect(annotype,"TE")) %>% select(-contains("dis"))

All_rnaplex3 <- bind_rows(gene_rnaplex3 %>% anti_join(gene_remove, by = "readname") %>% mutate(pairtype = "gene"),
                          TE_rnaplex3 %>% mutate(pairtype = "TE")) 

All_rnaplex3 %>% 
  write_tsv(paste0(outdir,"vienna/","All_model_vienna.tsv"))

All_rnaplex3 %>% 
  select(targetname,sstart2,send2) %>% 
  arrange(targetname,sstart2,send2) %>% 
  write_tsv(paste0(outdir,"vienna/","All_model_vienna_pre.bed"),col_names =  FALSE)

