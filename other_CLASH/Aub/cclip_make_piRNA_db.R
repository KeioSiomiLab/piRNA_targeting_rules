#######################################################################################################
# make piRNA database
#######################################################################################################

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
inputdir <- "/media/hermione/cCLIP/"
out <- "/media/hermione/cCLIP/"

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
  #file = "SRR8791716"
  #threshold =10
  #read file
  tx1 <- read_tsv(paste0(inputdir,file, ".fa"), col_names = "seq") %>%
    filter(!(str_detect(seq, "length"))) %>%
    group_by(seq) %>% summarise(count = n()) %>%
    mutate(CPM = count*1000000/sum(count))
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
fasta_make("SRR3051366", 0.5)
fasta_make("SRR3051369", 0.5)


fasta_make2 <- function(file){
  tx3 <- read_tsv(paste0(out,"piRNA_", file, ".tsv")) %>%
    mutate(type = file)
  return(tx3)
}
piRNA1 <- fasta_make2("SRR3051366")
piRNA2 <- fasta_make2("SRR3051369")

piRNAseq_ano <- tibble(type = c("SRR3051366", "SRR3051369"),
                       type2 = c("piRNA#1", "piRNA#2"))
piRNAall <- bind_rows(piRNA1, piRNA2) %>%
  left_join(piRNAseq_ano, by = "type")

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

