# database make
# /usr/bin/Rscript --silent --slave --vanilla /media/pericles/CLASH/database/get_gene_database_annotation.R;
#######################################################################################################
# annotation
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/CLASH/database/anno/")


suppressMessages(suppressWarnings(require(tidyverse)))
gene_repeat <- read_tsv("/media/pericles/CLASH/database/anno/repeat_gene.txt",
                        col_names = c("chr1", "start1", "end1", "gene_id1", "score1", "strand1", 
                                      "chr2", "start2", "end2", "gene_id2", "score2", "strand2")) %>% 
  select(-contains("score"))
gene_repeat2 <- gene_repeat %>% 
  mutate(start = if_else(strand1 == "+",start2-start1,end1-end2),
         end = if_else(strand1 == "+",end2 - start1,end1-start2),
         type = if_else(strand1==strand2, "TE_sense", "TE_anti")) %>% 
  mutate(start = if_else(start<0,0,start)) %>% 
  rename(name = gene_id1,transcript = gene_id2) %>% 
  select(name, transcript, start, end,type) %>% 
  filter(!(name %in% c("FBgn0267704","FBgn0287603")))

GTF <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_use.tsv")
gene <- GTF %>% 
  filter(type =="gene") %>% 
  mutate(gene_start = start,gene_end = end) %>%
  select(name, gene_start,gene_end)
exon <- GTF %>% 
  filter(type =="exon") %>% filter(name != "FBgn0002781") %>% 
  left_join(gene, by = "name")
exon2 <- exon %>% 
  mutate(new_start = if_else(direction=="+", start - gene_start, gene_end-end),
         new_end = if_else(direction=="+", end - gene_start+1,gene_end-start+1)) %>% 
  select(name,transcript, new_start,new_end)
intron <- exon2 %>% 
  group_by(name,transcript) %>% nest_legacy() %>% 
  mutate(start = map(data, function(x) {start = head(x$new_end, -1)}),
         end = map(data, function(x) {end = tail(x$new_start, -1)})) %>% 
  select(-data) %>% unnest_legacy() %>% 
  mutate(type = "intron") %>% 
  filter(!(name %in% c("FBgn0267704","FBgn0287603")))
CDS53 <- GTF %>% 
  filter(type %in% c("CDS","3UTR","5UTR")) %>% filter(name != "FBgn0002781") %>% 
  left_join(gene, by = "name") %>% 
  mutate(new_start = if_else(direction=="+", start - gene_start, gene_end-end),
         new_end = if_else(direction=="+", end - gene_start+1,gene_end-start+1)) %>% 
  select(name,transcript, new_start,new_end, type) %>% 
  rename(start = new_start,end = new_end) %>% 
  filter(!(name %in% c("FBgn0267704","FBgn0287603")))
flamhap1_anno <- read_tsv("/media/pericles/TEfind/flam/rmblast/flam_hap1_TE.bed",
                      col_names = c("name","start","end","transcript","score","strand")) %>% 
  mutate(name = "FBgn9999996", transcript = transcript %>% str_remove("_\\d+_\\d+"),
         type = if_else(strand == "+","TE_sense", "TE_anti")) %>% select(-score,-strand)
flamhap2_anno <- read_tsv("/media/pericles/TEfind/flam/rmblast/flam_hap2_TE.bed",
                      col_names = c("name","start","end","transcript","score","strand")) %>% 
  mutate(name = "FBgn9999997", transcript = transcript %>% str_remove("_\\d+_\\d+"),
         type = if_else(strand == "+","TE_sense", "TE_anti")) %>% select(-score,-strand)
l20Ahap1_anno <- read_tsv("/media/pericles/TEfind/flam/rmblast/20A_hap1_TE.bed",
                      col_names = c("name","start","end","transcript","score","strand")) %>% 
  mutate(name = "FBgn9999998", transcript = transcript %>% str_remove("_\\d+_\\d+"),
         type = if_else(strand == "+","TE_sense", "TE_anti")) %>% select(-score,-strand)
l20Ahap2_anno <- read_tsv("/media/pericles/TEfind/flam/rmblast/20A_hap2_TE.bed",
                      col_names = c("name","start","end","transcript","score","strand")) %>% 
  mutate(name = "FBgn9999999", transcript = transcript %>% str_remove("_\\d+_\\d+"),
         type = if_else(strand == "+","TE_sense", "TE_anti")) %>% select(-score,-strand)
finalanno <- bind_rows(CDS53,intron,gene_repeat2,flamhap1_anno,flamhap2_anno,l20Ahap1_anno,l20Ahap2_anno) %>% 
  select(name, start, end, transcript,type)
finalanno %>%
  write_tsv("/media/pericles/CLASH/database/anno/dm6_gene_ano.bed",col_names = FALSE)

#######################################################################################################
# newTE
#######################################################################################################


gene_repeat_rm <- read_tsv("/media/pericles/CLASH/database/anno/repeat_gene2.txt",
                        col_names = c("chr1", "start1", "end1", "gene_id1", "score1", "strand1", 
                                      "chr2", "start2", "end2", "gene_id2", "score2", "strand2")) %>% 
  select(-contains("score"))
gene_repeat2_rm <- gene_repeat_rm %>% 
  mutate(start = start2-start1,
         end = end2 - start1,
         type = if_else(strand1==strand2, "TE_sense", "TE_anti")) %>% 
  mutate(start = if_else(start<0,0,start)) %>% 
  rename(name = gene_id1,transcript = gene_id2) %>% 
  select(name, transcript, start, end,type)

gene_repeat2_rm %>% select(name, start, end, transcript,type) %>% 
  write_tsv("/media/pericles/CLASH/database/anno/dm6_gene_ano2.bed",col_names = FALSE)

gene_anno <- read_tsv("/media/pericles/CLASH/database/anno/dm6_gene.tsv") %>% 
  select(ID, name3) %>% 
  mutate(name = ID %>% str_remove("ID="))

gene_repeat2_rm2 <- gene_repeat2_rm %>% 
  left_join(gene_anno, by = "name") %>% 
  select(name3, start, end, transcript,type)

gene_repeat2_rm2 %>% 
  write_tsv("/media/pericles/CLASH/database/anno/dm6_gene_ano2_pre.bed",col_names = FALSE)
