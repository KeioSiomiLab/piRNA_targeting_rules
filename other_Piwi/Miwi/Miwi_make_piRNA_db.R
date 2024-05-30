# make piRNA database
# setwd("/media/pericles/CLASH/piRNA/")
  .libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

  inputdir <- "/media/hermione/otherdata_compare/Miwi/"
out <- "/media/hermione/otherdata_compare/Miwi/"
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(scales)))
suppressMessages(suppressWarnings(require(ggupset)))
suppressMessages(suppressWarnings(require(UpSetR)))
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
options(digits=7)

file_name_list <- c("SRR21528491", "SRR21528492","SRR21528493", "SRR21528494","SRR21528495","SRR21528497","SRR21528498", 
                    "SRR21528499","SRR21528500", "SRR21528501","SRR21528502","SRR21528503")

#spike in reads
#UGCUAGUCUGUUAUCGACCUGACCUCAUAG
#UGCUAGUCUGUUCGAUACCUGACCUCAUAG
#UGCUAGUCUGUUGUCACGAAGACCUCAUAG
#UGCUAGUCUUAUCGACCUCCUCAUAG
#UGCUAGUCUUCGAUACCUCCUCAUAG
#UGCUAGUCUUGUCACGAACCUCAUAG
#UGCUAGUUAUCGACCUUCAUAG
#UGCUAGUUCGAUACCUUCAUAG
#UGCUAGUUGUCACGAAUCAUAG

spikein <-tribble(
  ~seq, ~key,
  "UGCUAGUCUGUUAUCGACCUGACCUCAUAG", "spikein",
  "UGCUAGUCUGUUCGAUACCUGACCUCAUAG", "spikein",
  "UGCUAGUCUGUUGUCACGAAGACCUCAUAG", "spikein",
  "UGCUAGUCUUAUCGACCUCCUCAUAG", "spikein",
  "UGCUAGUCUUCGAUACCUCCUCAUAG", "spikein",
  "UGCUAGUCUUGUCACGAACCUCAUAG", "spikein",
  "UGCUAGUUAUCGACCUUCAUAG", "spikein",
  "UGCUAGUUCGAUACCUUCAUAG", "spikein",
  "UGCUAGUUGUCACGAAUCAUAG", "spikein") %>%mutate(seq = seq %>% str_replace_all("U","T"))

spikein_nom <-tribble(
  ~sample, ~cell_num, ~spike_mol,
  "SRR21528491", 31400, 370,
  "SRR21528492", 68900, 4000,
  "SRR21528493", 47200, 3000,
  "SRR21528494", 89100, 4000,
  "SRR21528495", 48900, 3000,
  "SRR21528497", 60300, 4000,
  "SRR21528498", 63100, 4000,
  "SRR21528499", 99700, 4000,
  "SRR21528500", 112500, 4000,
  "SRR21528501", 116400, 4000,
  "SRR21528502", 137000, 4000,
  "SRR21528503", 112000, 4000)

avogadro_threshold <- 1 /(6.02214076 * 100000)
# 1 molecule in cell, atto mol

fasta_make <- function(file){
  tx3 <- read_tsv(paste0(inputdir,file, ".fa"), col_names = "seq") %>%
    filter(!(str_detect(seq, "length"))) %>%
    group_by(seq) %>% summarise(count = n()) %>% ungroup() %>% mutate(sample = file) %>% 
    left_join(spikein, by = "seq") %>% left_join(spikein_nom, by = "sample")
  spike_in_mol <- tx3 %>% filter(key == "spikein") %>% pull(count) %>% sum()
  tx4 <- tx3 %>% 
    mutate(abattomol = count * spike_mol/(spike_in_mol * cell_num)) %>% filter(is.na(key)) %>% select(sample,seq,count,abattomol)
  tx4 %>% write_tsv(paste0(out, "otherdata/piRNA_", file, "_for_graph.tsv"))
  tx5 <- tx4 %>% filter(abattomol > avogadro_threshold) %>% 
    mutate(annotation= paste0("piRNA", row_number()),
           index = row_number())
  tx5 %>% write_tsv(paste0(out,"piRNA_", file, ".tsv"))
  #to fasta
  fasta <- tx5 %>% mutate(annotation2 = paste0(">", annotation)) %>%
    select(-count, -abattomol,-sample) %>%
    gather(seq, annotation2, key = "type", value = "fasta") %>%
    arrange(index, type) %>% select(fasta)
  fasta %>% write_tsv(paste0(out,"piRNA_", file, ".fa"), col_names = FALSE)
}

for (i in seq_along(file_name_list)) {
  fasta_make(file_name_list[[i]])
}


fasta_graph_raw <- function(file){
  tx3 <- read_tsv(paste0(out,"piRNA_", file, ".tsv")) %>%
    mutate(type = file)
  return(tx3)
}
piRNA1 <- fasta_graph_raw(file_name_list[[1]])
piRNA2 <- fasta_graph_raw(file_name_list[[2]])
piRNA3 <- fasta_graph_raw(file_name_list[[3]])
piRNA4 <- fasta_graph_raw(file_name_list[[4]])
piRNA5 <- fasta_graph_raw(file_name_list[[5]])
piRNA6 <- fasta_graph_raw(file_name_list[[6]])
piRNA7 <- fasta_graph_raw(file_name_list[[7]])
piRNA8 <- fasta_graph_raw(file_name_list[[8]])
piRNA9 <- fasta_graph_raw(file_name_list[[9]])
piRNA10 <- fasta_graph_raw(file_name_list[[10]])
piRNA11 <- fasta_graph_raw(file_name_list[[11]])
piRNA12 <- fasta_graph_raw(file_name_list[[12]])

library(GGally)
cor_func <- function(data,mapping,method,symbol){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x,y, method = method,use = 'complete.obs')
  colFn <- colorRampPalette(c("brown1","white","dodgerblue"),interpolate = "spline")
  fill <- colFn(100)[findInterval(corr,seq(-1,1,length = 100))]
  ggally_text(label = paste(symbol,as.character(round(corr, 2))),
              mapping = aes(),
              xP = 0.5,yP = 0.5,color = "black", size = 2.5)+
    theme(panel.background = element_rect(fill = fill))
}
theme_set(theme_minimal())
ppi <- 300
All_count <- bind_rows(piRNA1,piRNA2,piRNA3, piRNA4,piRNA5,piRNA6,piRNA7,piRNA8,piRNA9,piRNA10,piRNA11,piRNA12) %>% select(seq,sample,abattomol) %>% 
  pivot_wider(names_from = sample, values_from = abattomol, values_fill =list(abattomol=0))
png(paste0(out,"figures/","cor_each_sample.png"), width = 10*ppi, height = 10*ppi, res = ppi)
pm <- ggpairs(data = All_count %>% select(-seq) %>%
                mutate_at(vars(contains("SRR")),list(~log10(.*1000000+1))),
              xlab = "log10(target count + 1) ",ylab = "log10(target count + 1)",
              lower = list(continuous = wrap("points", alpha = 0.7, size = 0.2)),
              upper = list(continuous = wrap(cor_func, method = "pearson", symbol = expression('p ='))))+
  theme(text = element_text(size = 6), panel.grid.minor = element_blank())
pm2 <- pm
for(i in 2:pm$nrow) {
  for(j in 1:(i-1)) {
    pm2[i,j] <- pm[i,j] +
      scale_x_continuous(limits = c(0, 3.45)) +
      scale_y_continuous(limits = c(0, 3.45))}}
pm2
dev.off()


piRNAall <- bind_rows(piRNA1,piRNA2,piRNA3, piRNA4,piRNA5,piRNA6,piRNA7,piRNA8,piRNA9,piRNA10,piRNA11,piRNA12)

piRNAcheck <- piRNAall %>% select(seq, type) %>%
  group_by(seq) %>%
  summarize(count = length(type),
            type = list(type)) %>% ungroup() %>% 
  filter(str_count(seq, "A")/str_length(seq) <= 0.9) %>% 
  filter(str_count(seq, "T")/str_length(seq) <= 0.9) %>% 
  filter(str_count(seq, "G")/str_length(seq) <= 0.9) %>% 
  filter(str_count(seq, "C")/str_length(seq) <= 0.9)

fasta_all <- piRNAcheck %>% filter(count >= 6) %>%
  select(seq) %>%
  mutate(annotation= paste0("piRNA", row_number()),
         index = row_number(),
         annotation2 = paste0(">", annotation, "_piRNA_target_piRNA")) %>%
  gather(seq, annotation2, key = "type", value = "fasta") %>%
  arrange(index, type) %>% select(fasta)
fasta_all %>% write_tsv(paste0(out,"piRNA_all.fasta"), col_names = FALSE)

fasta_tsv <- piRNAcheck %>% filter(count >= 6) %>%
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

files3 <- list(piRNA1,piRNA2,piRNA3,piRNA4,piRNA5,piRNA6,piRNA7,piRNA8,piRNA9,piRNA10,piRNA11,piRNA12)
files4 <- file_name_list
df2 <- vector("list", length(files3))
for (i in seq_along(files4)) {
  df2[[i]] <- piRNA_tsv %>% 
    left_join(read_tsv(paste0(out, "otherdata/piRNA_", files4[[i]], "_for_graph.tsv")) %>% 
                select(seq, abattomol), by = "seq")
  df2[[i]][["sample"]] <- files4[[i]]
}


merge_CPM <- df2 %>% bind_rows() %>%
  mutate(abattomol = if_else(is.na(abattomol), 0, abattomol)) %>%
  select(name, seq, abattomol, sample, pilen) %>% 
  pivot_wider(names_from = sample, values_from = abattomol,values_fill =list(abattomol=0))

merge_CPM2 <- merge_CPM %>% mutate(abattomol = (SRR21528491+SRR21528492+SRR21528493+SRR21528494+SRR21528495+SRR21528497+SRR21528498+
                                                SRR21528499+SRR21528500+SRR21528501+SRR21528502+SRR21528503)/12)


merge_CPM2 %>%
  write_tsv(paste0(out,"piRNA_abattomol.tsv"))
  
