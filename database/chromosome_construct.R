# chr process
# /usr/bin/Rscript --silent --slave --vanilla chromosome_construct.R;
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
makefasta <- function(data, fasname, seqr){
  seqr <- rlang::enquo(seqr)
  x <- rlang::enquo(fasname)
  res <- data %>% 
    mutate(annotation2 = paste0(">", !!x), index = row_number()) %>% 
    select(annotation2,!!seqr,index) %>% 
    pivot_longer(cols = -index, names_to = "type", values_to = "fasta") %>% 
    arrange(index,type) %>% select(fasta)
  return(res)
}
#first TE annotation
genomedir <- commandArgs(trailingOnly = TRUE)[1]
# genomedir <- "/media/pericles/TEfind/database/dm6_process.tsv"

dm6 <- read_tsv(genomedir, col_names = c("name", "seq"))
dm6_2 <- dm6 %>% separate(name, sep = "\\s", into = "name2") %>% 
  filter(name2 %in% c("2L","2R","3L","3R","4","X","Y")) %>% 
  select(name2, seq) 
dm6_2 %>%  write_tsv("database/dm6.tsv", col_names = FALSE) %>% 
  makefasta(fasname = name2, seqr = seq) %>% 
  write_tsv("database/dm6.fasta", col_names = FALSE) 

dm6_chr <- dm6 %>% separate(name, sep = "\\s", into = "name2") %>% 
  mutate(name2 = if_else(name2 %in% c("2L","2R","3L","3R","4","X","Y"), paste0("chr", name2), name2)) %>% 
  select(name2, seq) 
dm6_chr %>% filter(name2 %in% paste0("chr", c("2L","2R","3L","3R","4","X","Y"))) %>% 
  write_tsv("database/dm6_chr.tsv", col_names = FALSE) %>% 
  makefasta(fasname = name2, seqr = seq) %>% 
  write_tsv("database/dm6_chr.fasta", col_names = FALSE) 
dm6_chr %>% write_tsv("database/dm6_chr_full.tsv", col_names = FALSE) %>% 
  makefasta(fasname = name2, seqr = seq) %>% 
  write_tsv("database/dm6_chr_full.fasta", col_names = FALSE) 

#each genome
tes <- dm6_chr %>% 
  mutate(index = row_number(), name3 = name2) %>% 
  mutate(annotation2 = paste0(">", name2)) %>% 
  gather(seq, annotation2, key = "type", value = "fasta") %>% 
  arrange(index, type) %>% select(name3,fasta)

tx <- tes %>% group_by(name3) %>% nest_legacy()
tx2 <- map(1:length(tx$name3), function(x) data.frame(tx$data[[x]]))
names(tx2) <- tx$name3
data <- map(1:length(tx2), function(x) write_tsv(tx2[[x]], 
                                                 paste0("database/allchromosome/", names(tx2[x]), ".fasta"),
                                                 col_names = FALSE))
