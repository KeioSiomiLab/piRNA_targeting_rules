# database make
# /usr/bin/Rscript --silent --slave --vanilla filter_endtoend.R endtoend/${X}_map2.sam Process/${X}_del.tsv Process/${X}_comp2.fasta;
#######################################################################################################
# annotation
#######################################################################################################
# samdir <- "/media/pericles/CLASH/endtoend/rep1_map2.sam"
# readdir <- "/media/pericles/CLASH/Process/rep1_del.tsv"
# outdir <- "/media/pericles/CLASH/Process/rep1_comp2.fasta"
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))

samdir <- commandArgs(trailingOnly = TRUE)[1]
readdir <- commandArgs(trailingOnly = TRUE)[2]
outdir <- commandArgs(trailingOnly = TRUE)[3]

sam1 <- read_tsv(samdir, 
                 col_names = c("name", "flag", "target", "start", "qual", "cigar")) %>% 
  drop_na() %>% 
  filter(flag != 4 & flag !=16)  %>% 
  filter(!(str_detect(target,"lncRNA:flam")))%>% 
  filter(!(str_detect(target,"lncRNA:20A")))
#  filter(!(str_detect(cigar, "I|D|S")))
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
