# database make
# /usr/bin/Rscript  cclip_filter_endtoend.R endtoend/{/.}_map2.sam Process/{/.}_del.tsv Process/{/.}_comp2.fasta
#######################################################################################################
# annotation
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
# samdir <- "/media/hermione/cCLIP/endtoend/SRR3051373_map2.sam"
# readdir <- "/media/hermione/cCLIP/Process/SRR3051373_del.tsv"
# outdir <- "/media/hermione/cCLIP/Process/SRR3051373_comp2.fasta"

suppressMessages(suppressWarnings(require(tidyverse)))

samdir <- commandArgs(trailingOnly = TRUE)[1]
readdir <- commandArgs(trailingOnly = TRUE)[2]
outdir <- commandArgs(trailingOnly = TRUE)[3]

sam1 <- read_tsv(samdir, 
                 col_names = c("name", "flag", "target", "start", "qual", "cigar")) %>% 
  drop_na() %>% 
  filter(flag != 4 & flag !=16) %>% 
  filter(!(str_detect(cigar, "I|D|S")))
sam2 <- sam1 %>% mutate(safer = "safer") %>% select(name, safer) %>% 
  distinct()
  
read1 <- read_tsv(readdir,col_names = c("name", "seq"))

read2 <- read1 %>% anti_join(sam2, by = "name")

makefasta <- function(data){
  res <- data %>% 
    mutate(name = paste0(">",name),index = row_number()) %>% 
    select(name,seq,index) %>% 
    pivot_longer(cols = - index, names_to = "type", values_to = "fasta") %>% 
    arrange(index,type) %>% select(fasta)
  return(res)
}
fasta1 <- read2 %>% makefasta()
fasta1 %>% write_tsv(outdir, col_names = FALSE)
