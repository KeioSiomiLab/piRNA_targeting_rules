# make piRNA database
# setwd("/media/pericles/CLASH/piRNA/")
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

inputdir <- "/media/hermione/otherdata_compare/Aub_testis/"
out <- "/media/hermione/otherdata_compare/Aub_testis/"

filename1 <- "SRR12213370"
filename2 <- "SRR12213371"

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
fasta_make(filename2, 0.5)

fasta_graph_raw <- function(file){
  tx3 <- read_tsv(paste0(out, "otherdata/piRNA_", file, "_for_graph.tsv")) %>%
    mutate(type = file)
  return(tx3)
}
piRNA1 <- fasta_graph_raw(filename1)
piRNA2 <- fasta_graph_raw(filename2)


piRNAseq_ano <- tibble(type = c(filename1, filename2),
                       type2 = c("piRNA#1", "piRNA#2"))
piRNAgraph <- bind_rows(piRNA1, piRNA2) %>%
  left_join(piRNAseq_ano, by = "type")
threshold <- 0.5
piRNAindex <- piRNAgraph %>% group_by(type2) %>%
  summarise(threshold = threshold * sum(count/1000000)) %>% ungroup() %>%
  mutate(threshold = if_else(threshold >= 3, threshold, 3))
piRNAgraph %>% group_by(type2, count) %>%
  summarise(count2 = n()) %>% ungroup() %>%
  ggplot() +
  geom_point(aes(x = count, y = count2, color = type2), size = 1, alpha = 0.7, shape= 19) +
  geom_vline(data = piRNAindex, aes(xintercept = threshold), color = "red") +
  geom_text(data = piRNAindex, aes(x = threshold , y = 100000, label = paste0(" threshold = ",round(threshold, digits = 1))),
            hjust=-0.2, size = 3) +
  xlab("piRNA enrichment of each data") +
  ylab("count of same enrichment") +
  scale_x_log10(breaks = 10 ^ (0:5),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = 10 ^ (0:6),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  facet_wrap( ~ type2)+
  theme_minimal()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  guides(color = FALSE) +
  ggsave(paste0(out,"figures/piRNA_threshold.png") , width =7, height =5, dpi = 300)

fasta_make2 <- function(file){
  tx3 <- read_tsv(paste0(out,"piRNA_", file, ".tsv")) %>%
    mutate(type = file)
  return(tx3)
}
piRNA1 <- fasta_make2(filename1)
piRNA2 <- fasta_make2(filename2)

piRNAseq_ano <- tibble(type = c(filename1, filename2),
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

files3 <- list(piRNA1,piRNA2)
files4 <- c(filename1, filename2)
df2 <- vector("list", length(files3))
for (i in seq_along(files4)) {
  df2[[i]] <- piRNA_tsv %>% left_join(read_tsv(paste0(out, "otherdata/piRNA_", files4[[i]], "_for_graph.tsv")) %>% 
                select(seq, count,CPM), by = "seq")
  df2[[i]][["sample"]] <- files4[[i]]
}
  
merge_CPM <- df2 %>% bind_rows() %>%
  mutate(CPM = if_else(is.na(CPM), 0, CPM)) %>%
  select(name, seq, CPM, sample, pilen) %>%
  pivot_wider(names_from = sample, values_from = CPM)

merge_CPM2 <- merge_CPM %>% mutate(CPM_avg = ( SRR12213370 + SRR12213371 ) / 2)


merge_CPM2 %>%
  write_tsv(paste0(out,"piRNA_CPM.tsv"))
  
