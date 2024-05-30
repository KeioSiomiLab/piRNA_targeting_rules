# clash_tf_sep2
# /usr/bin/Rscript --silent --slave --vanilla clash_tf_sep2_splice.R splice/${X}_tf_splice_map.sam splice/${X}_comp_tf_splice.tsv piRNA/anno/piRNA_all2.tsv splice/${X}_chimera;

#######################################################################################################
# separate chimera file into piRNA and target sequence in splicing reads
#######################################################################################################

# hybdir <- "/media/pericles/CLASH/splice/rep1_tf_splice_map.sam"
# CLASHdir <- "/media/pericles/CLASH/splice/rep1_comp_tf_splice.tsv"
# piRNAdir <- "/media/pericles/CLASH/piRNA/anno/piRNA_all2.tsv"
# outname <- "/media/pericles/CLASH/splice/rep1_chimera"
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

suppressMessages(suppressWarnings(require(tidyverse)))

hybdir <- commandArgs(trailingOnly = TRUE)[1]
CLASHdir <- commandArgs(trailingOnly = TRUE)[2]
piRNAdir <- commandArgs(trailingOnly = TRUE)[3]
outname <- commandArgs(trailingOnly = TRUE)[4]

hyb1 <- read_tsv(hybdir, 
                 col_names = c("target", "flag", "target2", "start", "qual", "cigar")) %>% 
  drop_na() %>% 
  filter(flag != 4 & flag !=16) %>% 
  filter((str_detect(cigar, "N"))) 

piRNA <- read_tsv(piRNAdir, col_names = c("piRNA", "piRNAseq"))

CLASH <- read_tsv(CLASHdir, col_names = c("target", "targetseq"))

hyb_all_pre <- hyb1 %>% left_join(CLASH, by="target") %>% 
  separate(target, sep = "@",into = c("piRNA", "readname", "direction")) %>% 
  left_join(piRNA, by = "piRNA")
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
  select(readname,piRNA,pilen,piRNAseq,target2,bedstart,bedend,targetlen,targetseq,cigar,direction)

hyb_all %>% write_tsv(paste0(outname, ".tsv"))

hyball2 <- hyb_all %>% select(readname,piRNA,piRNAseq,targetseq) %>% 
  separate(piRNA,sep = "_",into = c("piRNA", "dis","piRNAtype")) %>% distinct()

hyball3 <- hyball2 %>% mutate(index = row_number(),piRNA = paste0(">",piRNA), readname = paste0(">",readname)) %>% 
  pivot_longer(cols = c("piRNA","piRNAseq","readname","targetseq"), names_to = "type", values_to = "fasta") %>% 
  select(fasta)
hyball3 %>% write_tsv(paste0(outname, "_rnaplex.fasta"), col_names = FALSE)
