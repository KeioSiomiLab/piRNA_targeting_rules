# make piRNA database
# setwd("/media/pericles/CLASH/piRNA/")
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

inputdir <- "/media/hermione/otherdata_compare/Aub_embryo/"
out <- "/media/hermione/otherdata_compare/Aub_embryo/"

filename1 <- "SRR3051363"

suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(scales)))
suppressMessages(suppressWarnings(require(ggupset)))
suppressMessages(suppressWarnings(require(UpSetR)))
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
options(digits=7)

fasta_make <- function(file, threshold){
  tx1 <- read_tsv(paste0(inputdir,file, ".fa"), col_names = "seq") %>%
    filter(!(str_detect(seq, "length"))) %>%
    group_by(seq) %>% summarise(count = n()) %>%
    mutate(CPM = count*1000000/sum(count))
  tx1 %>% write_tsv(paste0(out, "otherdata/piRNA_", file, "_for_graph.tsv"))
  #filter and index
  tx3 <- tx1 %>% filter(CPM >= threshold) %>% filter(count >= 3) %>%
    mutate(annotation= paste0("piRNA", row_number()),
           index = row_number())
  tx3 %>% write_tsv(paste0(out,"piRNA_", file, ".tsv"))
  #to fasta
  fasta <- tx3 %>% mutate(annotation2 = paste0(">", annotation)) %>%
    select(-count, -CPM) %>%
    gather(seq, annotation2, key = "type", value = "fasta") %>%
    arrange(index, type) %>% select(fasta)
  fasta %>% write_tsv(paste0(out,"piRNA_", file, ".fa"), col_names = FALSE)
}
fasta_make(filename1, 0.5)

fasta_graph_raw <- function(file){
  tx3 <- read_tsv(paste0(out, "otherdata/piRNA_", file, "_for_graph.tsv")) %>%
    mutate(type = file)
  return(tx3)
}
piRNA1 <- fasta_graph_raw(filename1)

fasta_make2 <- function(file){
  tx3 <- read_tsv(paste0(out,"piRNA_", file, ".tsv")) %>%
    mutate(type = file)
  return(tx3)
}
piRNA1 <- fasta_make2(filename1)

piRNAall <- piRNA1

piRNAcheck <- piRNAall %>% select(seq, type) %>%
  group_by(seq) %>%
  summarize(count = length(type),
            type = list(type))

fasta_all <- piRNAcheck %>% filter(count >= 1) %>%
  select(seq) %>%
  mutate(annotation= paste0("piRNA", row_number()),
         index = row_number(),
         annotation2 = paste0(">", annotation, "_piRNA_target_piRNA")) %>%
  gather(seq, annotation2, key = "type", value = "fasta") %>%
  arrange(index, type) %>% select(fasta)
fasta_all %>% write_tsv(paste0(out,"piRNA_all.fasta"), col_names = FALSE)

fasta_tsv <- piRNAcheck %>% filter(count >= 1) %>%
  select(seq) %>%
  mutate(annotation= paste0("piRNA", row_number()),
         index = row_number(),
         name = paste0( annotation, "_piRNA_target_piRNA")) %>% 
  select(name,seq)
fasta_tsv %>% write_tsv(paste0(out,"piRNA_all.tsv"), col_names = FALSE)
##########################################################
#CPM cal
##########################################################

piRNA_tsv <- fasta_tsv %>%
  mutate(name = name %>% str_extract("piRNA\\d+"),
         pilen = str_length(seq))

merge_CPM <- piRNA_tsv %>% left_join(piRNA1 %>% select(seq, count,CPM), by = "seq") %>%
  mutate(CPM = if_else(is.na(CPM), 0, CPM),
         sample = "SRR3051363") %>%
  select(name, seq, CPM, sample, pilen) %>%
  pivot_wider(names_from = sample, values_from = CPM)

merge_CPM2 <- merge_CPM %>% mutate(CPM_avg = SRR3051363)


merge_CPM2 %>%
  write_tsv(paste0(out,"piRNA_CPM.tsv"))
  
