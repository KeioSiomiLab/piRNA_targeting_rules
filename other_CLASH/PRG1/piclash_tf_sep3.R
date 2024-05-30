#######################################################################################################
# separate chimera file into piRNA and target sequence
#######################################################################################################
# /usr/bin/Rscript  piclash_tf_sep3.R Process3/{/.}_tf_map.sam Process2/{/.}_pimap.tsv Process3/{/.}_chimera
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)

# hybdir <- "/media/hermione/PRG1_CLASH/Process3/SRR6512652_tf_map.sam"
# CLASHdir <- "/media/hermione/PRG1_CLASH/Process2/SRR6512652_pimap.tsv"
# outname <- "/media/hermione/PRG1_CLASH/Process3/SRR6512652_chimera"

hybdir <- commandArgs(trailingOnly = TRUE)[1]
CLASHdir <- commandArgs(trailingOnly = TRUE)[2]
outname <- commandArgs(trailingOnly = TRUE)[3]

hyb1 <- read_tsv(hybdir,
                 col_names = c("name", "flag", "target2", "start", "qual", "cigar")) %>%
  drop_na() %>%
  filter(flag != 4)

CLASH <- read_tsv(CLASHdir)

hyb_all_pre <- hyb1 %>% left_join(CLASH, by="name") %>%
  rename(targetlen = left_length,
         readname = name,
         piRNAcount = count,
         targetseq = seq2) %>%
  mutate(pilen = str_length(piRNAseq),
         bedstart = start-1,
         bedend =start + str_length(targetseq)-1,
         piRNA = readname)
cigar_index <- hyb_all_pre %>% select(cigar) %>% distinct() %>%
  mutate(pack = map(cigar,function(x){str_split(x,pattern = "[:upper:]")}),
         titania = map(cigar,function(x){str_split(x,pattern = "[:digit:]")})) %>%
  unnest_legacy() %>%
  mutate(pack2 = map(pack, function(x) {x[which(x != "")]}),
         titania2 = map(titania, function(x) {x[which(x != "")]}),
         titania3 = map(titania2, function(x) {(which(x=="I"|x=="S"))}),
         end = map2(pack2, titania3, function(x,y) {replace(x,y,0) %>% as.integer() %>% sum()})) %>%
  select(-contains("pack"),-contains("titania")) %>%
  unnest_legacy()
hyb_all <- hyb_all_pre %>% left_join(cigar_index,by = "cigar") %>%
  mutate(pilen = str_length(piRNAseq),
         bedstart = start-1,
         bedend =start + end-1,
         target2 = str_remove_all(target2, "\\-"),
         targetlen = str_length(targetseq)) %>%
#  filter(targetlen > pilen) %>%
  select(readname,piRNA,pilen,piRNAseq,piRNAcount,target2,bedstart,bedend,targetlen,targetseq,cigar)

hyb_all %>% write_tsv(paste0(outname, ".tsv"))

hyball2 <- hyb_all %>% select(readname,piRNA,piRNAseq,targetseq) %>%
  distinct()

hyball3 <- hyball2 %>% mutate(index = row_number(),piRNA = paste0(">piRNA",piRNA), readname = paste0(">",readname)) %>%
  pivot_longer(cols = c("piRNA","piRNAseq","readname","targetseq"), names_to = "type", values_to = "fasta") %>%
  select(fasta)
hyball3 %>% write_tsv(paste0(outname, "_rnaplex.fasta"), col_names = FALSE)
