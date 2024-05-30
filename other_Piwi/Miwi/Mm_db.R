
############################################################################################################################################################################
# GTF process
############################################################################################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
options(scipen=100)
#gtf <- read_tsv("/media/hermione/otherdata_compare/Siwi/database/Bomo_gene_models.gtf",
#                skip = 5,
#                col_names= c("chr", "source", "type", "start", "end", "len", "direction", "score", "gene_name"))

gtf <- read_tsv("/media/qqprosperodd/miranda/Refee/database/gencode.vM25.annotation.gtf",
                skip = 5,
                col_names= c("chr", "source", "type", "start", "end", "len", "direction", "score", "gene_name"))

annotation1 <- gtf %>% filter(type =="gene") %>%
  select(gene_name) %>% #filter(str_detect(gene_name, "LOC\\d+")) %>% 
  mutate(gene_id = str_extract(gene_name,'gene_id "[:graph:]+"') %>% str_remove("gene_id ") %>% str_remove_all('"'),
         gene_type = str_extract(gene_name,'gene_type "[:graph:]+"') %>% str_remove("gene_type ") %>% str_remove_all('"'),
         gene_symbol = str_extract(gene_name,'gene_name "[:graph:]+"') %>% str_remove("gene_name ") %>% str_remove_all('"'))

annotation1 %>% 
  write_tsv("/media/hermione/otherdata_compare/Miwi/database/Mm_gene_annotation.txt")

#

gtf %>%  filter(type =="gene") %>%
  mutate(gene_id = str_extract(gene_name,'gene_id "[:graph:]+"') %>% str_remove("gene_id ") %>% str_remove_all('"')) %>% 
  mutate(strand = direction,
                start = start -1,
                score = 100) %>% 
  select(chr, start, end, gene_id, score, strand) %>% 
  arrange(chr,start,end) %>% 
  write_tsv("/media/hermione/otherdata_compare/Miwi/database/gene.bed",col_names = FALSE)

finalanno <- gtf %>% filter(type == "exon") %>% 
  mutate(gene_id = str_extract(gene_name,'gene_id "[:graph:]+"') %>% str_remove("gene_id ") %>% str_remove_all('"'),
         gene_type = str_extract(gene_name,'gene_type "[:graph:]+"') %>% str_remove("gene_type ") %>% str_remove_all('"')) %>% 
  select(chr, start, end, gene_id,gene_type)
finalanno %>%
  write_tsv("/media/hermione/otherdata_compare/Miwi/database/Mm_gene_ano.bed",col_names = FALSE)


#repeatmasker
filedir <- "/media/hermione/otherdata_compare/Miwi/database/RepMas/mm10.fa.out.gff"
pbsv_repmask <- read_tsv(filedir,  col_names = c("name2","source", "dis1","insstart","insend", "score","strand","dis2","name"),skip = 3)

pbsv_repmask2 <- pbsv_repmask %>%
  separate(name,into = c("dis3","TEname","TEstart","TEend"), sep = " ") %>%
  mutate(TEname = TEname %>% str_remove_all('"|Motif:'))

pbsv_repmask3 <- pbsv_repmask2 %>%
  mutate(TEstrand = if_else(strand == "+", "TE_sense","TE_anti")) %>% 
  select(name2,insstart,insend,TEname,TEstrand,strand)%>% distinct() 
pbsv_repmask3 %>%
  write_tsv("/media/hermione/otherdata_compare/Miwi/database/Mm_TE_ano.bed",col_names = FALSE)



