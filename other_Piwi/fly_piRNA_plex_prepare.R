

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))

outdir <- "/media/hermione/otherdata_compare/Aub_testis/"

outdir <- commandArgs(trailingOnly = TRUE)[1]

#######################################################################################################
# processing piRNA mapping to the models
#######################################################################################################

piRNA_seq <- read_tsv(paste0(outdir, "piRNA_CPM.tsv")) %>% 
  select(name,seq, pilen,CPM_avg)

piRNA_model <- read_tsv(paste0(outdir,"blastout/","piRNA.txt"),
                       col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                     "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>% 
  mutate(qseqid = qseqid %>% str_extract("piRNA\\d+")) %>%
  left_join(piRNA_seq %>% rename(qseqid = name), by = "qseqid") %>% 
  rename(CPM = CPM_avg) %>% 
  mutate(type = case_when((pilen == length & gapopen == 0) & mismatch==0 ~ "perfect match",
                          (qstart %in% c(1,2) & qend >= 20) & (gapopen == 0 & mismatch==0) ~ "functional",
                          TRUE ~ "not functional"))

TE_seq <- read_tsv("/media/pericles/CLASH/database/reffasta/ensembl.tsv", col_names = c("name", "seq")) 

piRNA_model2 <-  piRNA_model %>% 
 # filter(type != "not functional") %>% 
  mutate(strand = if_else(sstart > send, "antisense", "sense"),
         sstart2 = if_else(strand == "antisense", send, sstart),
         send2 = if_else(strand == "antisense", sstart, send)) %>% 
  left_join(TE_seq %>% rename(sseqid = name, targetseq = seq), by = "sseqid") %>% 
  mutate(targetseq = targetseq %>% str_sub(start = sstart2-5, end = sstart2 + 35),
         targetlen = str_length(targetseq)) %>% 
  filter(strand == "antisense") %>% filter(targetlen != 0) %>% 
  mutate(readname = paste0("piRNAmodel",row_number())) %>% 
  rename(piRNAname = qseqid, targetname = sseqid, piRNAseq = seq) %>% 
  select(readname, piRNAname, targetname, pilen,piRNAseq,qstart,qend,targetlen,targetseq,sstart2,send2,evalue,bitscore, CPM,type, strand)
  
piRNA_model2 %>% write_tsv(paste0(outdir,"rnaplex/", "model_RNAplex.tsv"))

piRNA_model3 <- piRNA_model2 %>% mutate(index = row_number(),piRNA = paste0(">",piRNAname), readname = paste0(">",readname)) %>% 
  select(index, piRNA, readname, piRNAseq,targetseq) %>% 
  pivot_longer(cols = c("piRNA","piRNAseq","readname","targetseq"), names_to = "type", values_to = "fasta") %>% 
  select(fasta)
piRNA_model3 %>% write_tsv(paste0(outdir,"rnaplex/", "model_RNAplex.fasta"), col_names = FALSE)

#######################################################################################################
# processing piRNA mapping to the models
#######################################################################################################
# TE mapped piRNAs
# remove TE mapped region
remove_piRNA <- piRNA_model2 %>% 
  select(piRNAname,type) %>% rename(remove_type = type) %>% distinct()

piRNA_model_gene <- read_tsv(paste0(outdir,"blastout/","gene.txt"),
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                      "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>% 
  mutate(qseqid = qseqid %>% str_extract("piRNA\\d+")) %>%
  left_join(piRNA_seq %>% rename(qseqid = name), by = "qseqid") %>% 
  rename(CPM = CPM_avg) %>% 
  mutate(type = case_when((pilen == length & gapopen == 0) & mismatch==0 ~ "perfect match",
                          (qstart %in% c(1,2) & qend >= 20) & (gapopen == 0 & mismatch==0) ~ "functional",
                          TRUE ~ "not functional"))
gene_seq <- read_tsv("/media/hermione/piRNAmodel/gene.tsv", col_names = c("name", "seq")) %>% 
  mutate(name = name %>% str_extract("FBgn\\d+"))


piRNA_model2_gene <-  piRNA_model_gene %>% 
 # filter(type != " not functional") %>% 
  left_join(remove_piRNA %>% rename(qseqid = piRNAname), by = "qseqid") %>% filter(is.na(remove_type)) %>% 
  mutate(strand = if_else(sstart > send, "antisense", "sense"),
         sstart2 = if_else(strand == "antisense", send, sstart),
         send2 = if_else(strand == "antisense", sstart, send)) %>% 
  left_join(gene_seq %>% rename(sseqid = name, targetseq = seq), by = "sseqid") %>% 
  mutate(targetseq = targetseq %>% str_sub(start = sstart2-5, end = sstart2 + 35),
         targetlen = str_length(targetseq)) %>% 
  filter(strand == "antisense") %>% filter(targetlen != 0) %>% 
  mutate(readname = paste0("piRNAmodel",row_number())) %>% 
  rename(piRNAname = qseqid, targetname = sseqid, piRNAseq = seq) %>% 
  select(readname, piRNAname, targetname, pilen,piRNAseq,qstart,qend,targetlen,targetseq,sstart2,send2,evalue,bitscore, CPM,type, strand)

piRNA_model2_gene %>% write_tsv(paste0(outdir,"rnaplex/", "model_RNAplex2.tsv"))

piRNA_model3_gene <- piRNA_model2_gene %>% mutate(index = row_number(),piRNA = paste0(">",piRNAname), readname = paste0(">",readname)) %>% 
  select(index, piRNA, readname, piRNAseq,targetseq) %>% 
  pivot_longer(cols = c("piRNA","piRNAseq","readname","targetseq"), names_to = "type", values_to = "fasta") %>% 
  select(fasta)
piRNA_model3_gene %>% write_tsv(paste0(outdir,"rnaplex/", "model_RNAplex2.fasta"), col_names = FALSE)


remove_gene_region <- piRNA_model2_gene %>% 
  mutate(start3 = sstart2,end3 = sstart2 + pilen) %>% 
  select(targetname, start3,end3, readname) %>% 
  write_tsv(paste0(outdir,"rnaplex/", "gene_remove1.bed"), col_names = FALSE)
  
