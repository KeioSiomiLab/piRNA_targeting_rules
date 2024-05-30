

#######################################################################################################
# chimera processing 2
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
# /usr/bin/Rscript chimera_process2.R overlap_res/{/.}_gene_overlap.bed overlap_res/{/.}_TE_overlap.bed chimera_res/{/.}_chimera.tsv


suppressMessages(suppressWarnings(require(tidyverse)))

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

# resfile <- "/media/hermione/CAGE/overlap_res/siEGFPrep1_1_gene_overlap.bed"
# junctionfile  <- "/media/hermione/CAGE/overlap/siEGFPrep1_1_junction.result"
# outfile <- "/media/hermione/CAGE/chimera_res/siEGFPrep1_1_chimera.tsv"


resfile <- commandArgs(trailingOnly = TRUE)[1]
junctionfile <- commandArgs(trailingOnly = TRUE)[2]
outfile <- commandArgs(trailingOnly = TRUE)[3]

res <- read_tsv(resfile,
                col_names = FALSE) %>% mutate(type = "gene")

junction_res <- read_tsv(junctionfile)

annotation <- read_tsv("/media/pericles/CLASH/database/anno/dm6_gene_annotation.tsv")
result <- junction_res %>% 
  left_join(res %>% select(X4, X10,X6) %>% rename(gene_name = X4,index=X10,strand = X6), by = "index") %>% 
  left_join(annotation %>% rename(gene_name = name3), by = "gene_name") %>% 
  mutate(symbol = if_else(is.na(symbol), "notannotated", symbol),
         rnatype = case_when(str_detect(symbol, "snRNA|lncRNA|rRNA|snoRNA|tRNA") ~ "rmstRNA",
                             TRUE ~ "normal"))

result2 <- result %>% 
  group_by(index) %>% filter(n() == 2)

rnatypeorder <- c("rmstRNA","normal")
miranda <- result %>% select(index,symbol,rnatype) %>%
  distinct() %>% mutate(rnatype = rnatype %>% factor(levels= rnatypeorder)) %>% 
  group_by(index) %>% slice_max(n=1, desc(rnatype), with_ties = TRUE) %>% 
  slice_max(n=1, desc(symbol), with_ties = FALSE) %>% ungroup() %>%
  select(index,symbol) 

chimaeric_junction1 <- miranda %>% left_join(result, by = c("index","symbol"))

chimaeric_junction1 %>% write_tsv(outfile)
