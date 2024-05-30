

#######################################################################################################
# mdg1 simulation piRNA selection
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
piRNA_seq <- read_tsv("/media/pericles/CLASH/piRNA/piRNA_CPM.tsv") %>%
  select(name,seq, pilen,CPM_avg)

piRNA_mdg1 <- read_tsv("/media/hermione/piRNA_luc/TE/mdg1.txt",
                       col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                     "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  mutate(qseqid = qseqid %>% str_extract("piRNA\\d+")) %>%
  left_join(piRNA_seq %>% rename(qseqid = name), by = "qseqid") %>% filter(sstart < send)


outdir <- "/media/hermione/piRNA_luc/figure/"


# 0.6% of total piRNA

anti_vector <- read_tsv("/media/hermione/piRNA_luc/mdg1_anti_mut.tsv",
                        col_names = c("vector","seq2"))
anti_vector2 <- anti_vector %>% mutate(new = list(piRNA_mdg1 %>% select(qseqid, seq) %>% distinct())) %>%
  unnest_legacy()



anti_vector4 <- anti_vector2 %>% mutate(index = row_number(),vector2 = vector, vector = paste0(">",vector), qseqid = paste0(">",qseqid)) %>%
  pivot_longer(cols = c("qseqid","seq","vector","seq2"), names_to = "type", values_to = "fasta")

tes <- anti_vector4 %>% select(vector2,fasta)

tx <- tes %>% group_by(vector2) %>% nest_legacy()
tx2 <- map(1:length(tx$vector2), function(x) data.frame(tx$data[[x]]))
names(tx2) <- tx$vector2
data <- map(1:length(tx2), function(x) write_tsv(tx2[[x]],
                                                 paste0("/media/hermione/piRNA_luc/","rnaplex_temp/", names(tx2[x]), ".fasta"),
                                                 col_names = FALSE))
anti_vector2 %>% mutate(vector2 = vector) %>%
  group_by(vector2) %>%
  group_walk(~ write_tsv(.x, paste0("/media/hermione/piRNA_luc/hyb_temp/",.y$vector2,".tsv")))

#######################################################################################################
# mdg1 and others
#######################################################################################################
piRNA_seq <- read_tsv("/media/pericles/CLASH/piRNA/piRNA_CPM.tsv") %>%
  select(name,seq, pilen,CPM_avg)

#plasmid
piRNA_plasmid <- read_tsv("/media/hermione/piRNA_luc/plasmid/plasmid.txt",
                       col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                     "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  mutate(qseqid = qseqid %>% str_extract("piRNA\\d+")) %>%
  left_join(piRNA_seq %>% rename(qseqid = name), by = "qseqid") %>% filter(sstart > send)

anti_vector_plasmid <- read_tsv("/media/hermione/piRNA_luc/plasmid/piRNA_plasmid_target.tsv",
                        col_names = c("vector","seq2"))
anti_vector2_plasmid <- anti_vector_plasmid %>% mutate(new = list(piRNA_plasmid %>% select(sseqid,qseqid, seq) %>% distinct())) %>%
  unnest_legacy() %>% filter(sseqid == vector)  %>% mutate(vector= vector %>% str_replace("piRNA","pi_RNA")) %>% 
  select(vector, seq2,qseqid,seq)

anti_vector4_plasmid <- anti_vector2_plasmid %>%
  mutate(index = row_number(),vector2 = vector, vector = paste0(">",vector), qseqid = paste0(">",qseqid)) %>%
  pivot_longer(cols = c("qseqid","seq","vector","seq2"), names_to = "type", values_to = "fasta")

tes_plasmid <- anti_vector4_plasmid %>% select(vector2,fasta)

tx_plasmid <- tes_plasmid %>% group_by(vector2) %>% nest_legacy()
tx2_plasmid <- map(1:length(tx_plasmid$vector2), function(x) data.frame(tx_plasmid$data[[x]]))
names(tx2_plasmid) <- tx_plasmid$vector2
data_plasmid <- map(1:length(tx2_plasmid), function(x) write_tsv(tx2_plasmid[[x]],
                                                 paste0("/media/hermione/piRNA_luc/","rnaplex_temp2/", names(tx2_plasmid[x]), ".fasta"),
                                                 col_names = FALSE))
anti_vector2_plasmid %>% mutate(vector2 = vector) %>%
  group_by(vector2) %>%
  group_walk(~ write_tsv(.x, paste0("/media/hermione/piRNA_luc/hyb_temp2/",.y$vector2,".tsv")))

#flam
piRNA_flam <- read_tsv("/media/hermione/piRNA_luc/flam/flam_anti.txt",
                          col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  mutate(qseqid = qseqid %>% str_extract("piRNA\\d+")) %>%
  left_join(piRNA_seq %>% rename(qseqid = name), by = "qseqid") %>% mutate(sseqid = sseqid %>% str_replace("\\.","_")) %>% filter(sstart > send)

anti_vector_flam <- read_tsv("/media/hermione/piRNA_luc/plasmid/flam_anti.tsv",
                                col_names = c("vector","seq2")) %>% mutate(vector = vector %>% str_replace("\\.","_"))
anti_vector2_flam <- anti_vector_flam %>% mutate(new = list(piRNA_flam %>% filter(send < sstart) %>% 
     mutate(send2 = send - (pilen-qend) - 5, sstart2 = sstart + qstart -1 + 5) %>% 
    select(sseqid,qseqid, seq,qstart,qend,sstart2,send2) %>% distinct())) %>%
  unnest_legacy() %>% filter(sseqid == vector) %>% 
  mutate(totallen = str_length(seq2) %>% as.double(),send2 = if_else(send2 < 0, 0, send2),
         sstart2 = if_else(sstart2 > totallen,totallen,sstart2)) %>% 
  mutate(seq2 = str_sub(seq2,start = send2, end = sstart2)) %>% select(vector, seq2,qseqid,seq)

anti_vector4_flam <- anti_vector2_flam %>% mutate(index = row_number(),vector2 = vector, vector = paste0(">",vector), qseqid = paste0(">",qseqid)) %>%
  pivot_longer(cols = c("qseqid","seq","vector","seq2"), names_to = "type", values_to = "fasta")

tes_flam <- anti_vector4_flam %>% select(vector2,fasta)

tx_flam <- tes_flam %>% group_by(vector2) %>% nest_legacy()
tx2_flam <- map(1:length(tx_flam$vector2), function(x) data.frame(tx_flam$data[[x]]))
names(tx2_flam) <- tx_flam$vector2
data_flam <- map(1:length(tx2_flam), function(x) write_tsv(tx2_flam[[x]],
                                                 paste0("/media/hermione/piRNA_luc/","rnaplex_temp2/", names(tx2_flam[x]), ".fasta"),
                                                 col_names = FALSE))
anti_vector2_flam %>% mutate(vector2 = vector) %>%
  group_by(vector2) %>%
  group_walk(~ write_tsv(.x, paste0("/media/hermione/piRNA_luc/hyb_temp2/",.y$vector2,".tsv")))









