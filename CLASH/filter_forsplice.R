# filter mapped reads from fasta
# /usr/bin/Rscript --silent --slave --vanilla filter_forsplice.R Process3/${X}_chimera.tsv Process3/${X}_comp_tf.tsv splice/${X}_comp_tf_splice.tsv;
#######################################################################################################
# annotation
#######################################################################################################
# setwd("/media/pericles/CLASH/")
# samdir <- "/media/pericles/CLASH/Process3/rep1_chimera.tsv"
# readdir <- "/media/pericles/CLASH/Process3/rep1_comp_tf.tsv"
# outdir <- "/media/pericles/CLASH/splice/rep1_comp_tf_splice.tsv"
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

suppressMessages(suppressWarnings(require(tidyverse)))

samdir <- commandArgs(trailingOnly = TRUE)[1]
readdir <- commandArgs(trailingOnly = TRUE)[2]
outdir <- commandArgs(trailingOnly = TRUE)[3]

sam1 <- read_tsv(samdir) %>% mutate(safer = "safer") %>%
  select(readname,safer) %>% distinct() %>% rename(name = readname)

read1 <- read_tsv(readdir,col_names = c("name", "seq")) %>% 
  separate(name, into = c("piRNA", "name", "direction"),sep = "@")

read2 <- read1 %>% anti_join(sam1, by = "name") %>% 
  unite(piRNA,name,direction,sep = "@",col="name")

read2 %>% write_tsv(outdir, col_names = FALSE)
