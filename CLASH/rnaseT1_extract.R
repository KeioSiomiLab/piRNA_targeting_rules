# extract RNAseT1 motif ... for memory usage
# /usr/bin/Rscript --silent --slave --vanilla rnaseT1_extract.R vienna/rep2

#######################################################################################################
# extract RNAseT1 motif
#######################################################################################################
# file_vienna <- "/media/pericles/CLASH/vienna/rep2"
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

suppressMessages(suppressWarnings(require(tidyverse)))

file_vienna <- commandArgs(trailingOnly = TRUE)[1]

gene_all <- read_tsv("/media/pericles/CLASH/database/dm6_gene.tsv",
                     col_names = c("target", "seq")) %>%
  separate(target, into = c("target","dis6","dis7"), sep = "_") %>% select(-dis6,-dis7) %>%
  mutate(target = target %>% str_remove("\\-"))

vienna <- read_tsv(paste0(file_vienna,"_vienna.tsv")) 

sub_RNAseT1_seq <- vienna %>% 
  select(target,bedstart) %>% distinct() %>% 
  left_join(gene_all,by = "target") %>% 
  mutate(seq = str_sub(seq,start = bedstart,end=bedstart) %>% str_to_upper())

vienna2 <- vienna %>% 
  left_join(sub_RNAseT1_seq,by = c("target","bedstart")) 

vienna2 %>% write_tsv(paste0(file_vienna, "_vienna2.tsv"))

