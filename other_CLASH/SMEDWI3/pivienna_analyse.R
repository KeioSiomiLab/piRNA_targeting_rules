#######################################################################################################
# make table
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

setwd("/media/hermione/SMEDWI3/vienna/")
library(tidyverse)
library(gplots)
library(RColorBrewer)
library(gtools)
library(gridExtra)
library(stringi)
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
filedir <- "/media/hermione/SMEDWI3/vienna/"
outdir <- "/media/hermione/SMEDWI3/"

############################################################
#set factor
#this is used for select read-target.
targetRNAlevel <- c("rRNA", "tRNA", "snRNA", "snoRNA", "mtrRNA", "miscRNA", "Ig", "prtr","mRNA",
                    "pseudo", "lincRNA",  "microRNA")

files3 <- dir(filedir, pattern = "_vienna.tsv$", full.names = TRUE)
files4 <- files3 %>%
  str_replace_all("_vienna.tsv$", "") %>%
  str_replace(paste0(filedir, "/"), "")
df2 <- vector("list", length(files3))
for (i in seq_along(files3)) {
  df2[[i]] <- read_tsv(files3[[i]], col_types = cols(dG = col_double()))
  df2[[i]][["sample"]] <- files4[[i]]
}
Vienna1 <- df2 %>% bind_rows() %>%
  mutate(dG = if_else(!(str_detect(piRNAmatch, "\\(")), Inf, dG))%>% 
  filter(str_count(targetseq, "A")/str_length(targetseq) <= 0.9) %>% 
  filter(str_count(targetseq, "T")/str_length(targetseq) <= 0.9) %>% 
  filter(str_count(targetseq, "G")/str_length(targetseq) <= 0.9) %>% 
  filter(str_count(targetseq, "C")/str_length(targetseq) <= 0.9) %>% 
  filter(str_count(piRNAseq, "A")/str_length(piRNAseq) <= 0.9) %>% 
  filter(str_count(piRNAseq, "T")/str_length(piRNAseq) <= 0.9) %>% 
  filter(str_count(piRNAseq, "G")/str_length(piRNAseq) <= 0.9) %>% 
  filter(str_count(piRNAseq, "C")/str_length(piRNAseq) <= 0.9) 

Vienna1 %>% write_csv(paste0(outdir, "vienna1fig/clash.csv"))

Vienna_select <- Vienna1 %>%
  write_csv(paste0(outdir, "vienna1fig/clash_selection.csv"))

Vienna_sle <- Vienna1 %>%
  select(sample, readname,dG, piRNA,pilen, piRNAseq,
         targetlen,targetseq,piRNAmatch, targetmatch) %>% distinct() %>% 
  group_by(sample, readname) %>%
  nest_legacy() %>%
  mutate(data2 = map(data, function(x) {x %>% filter(dG == min(dG)) %>%
      filter(abs(pilen-31.9)==min(abs(pilen-31.9)))})) %>%
  select(-data) %>% unnest_legacy()
Vienna_sle %>% write_csv(paste0(outdir, "vienna1fig/clash_colapse_pre.csv"))

#########################

Vienna_sle <- read_csv(paste0(outdir, "vienna1fig/clash_colapse_pre.csv"), col_types = cols(dG = col_double()))

#rRNA
rRNA <- read_tsv(paste0(outdir, "rRNA_filter/rRNA_piRNA.sam"), 
                 col_names = FALSE)

piRNA_fil <- left_join(Vienna_sle,rRNA %>% filter(X2==0) %>% select(X1,X2) %>% rename(piRNA=X1), by = "piRNA")

rRNA2 <- read_tsv(paste0(outdir, "database/raw/piRNA_rRNA.txt"), 
                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
rRNA3 <- read_tsv(paste0(outdir, "database/raw/piRNA_rRNA_silva.txt"), 
                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
# genome_unmapped piRNA
piRNA_un_genome <- read_tsv(paste0(outdir, "rRNA_filter/piRNA_genome_remove.sam"), 
                            col_names = FALSE)

Vienna_rRNA_remove <- Vienna_sle %>% 
  filter(targetlen > pilen) %>% 
  anti_join(rRNA2 %>% select(qseqid) %>% distinct() %>% rename(piRNA = qseqid), by = "piRNA") %>% 
  anti_join(rRNA3 %>% select(qseqid) %>% distinct() %>% rename(piRNA = qseqid), by = "piRNA") %>% 
  anti_join(piRNA_un_genome %>% select(X1) %>% distinct() %>% rename(piRNA = X1), by = "piRNA")

Vienna_rRNA_remove %>% write_csv(paste0(outdir, "vienna1fig/clash_colapse.csv"))
