# VCF file analyse
# remove repeat region that doesnt exist in OSC genome.
# /usr/bin/Rscript --silent --slave --vanilla vcf_process.R;
#######################################################################################################
# annotation
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
#######################################################################################################
# pbsv
#######################################################################################################
pbsv_vcf <- read_tsv("/media/pericles/TEfind/svcall/pacbio_pbsv.vcf",
                     skip = 34)
pbsv_vcf_lonins <- pbsv_vcf %>% filter(str_detect(ID, "INS")) %>%
  mutate(len = str_length(ALT)) %>%
  filter(len >50)
pbsv_vcf_lonins2 <- pbsv_vcf_lonins %>%
  mutate(start = POS, end = POS, strand = "+",name = paste0(ID,"@", FILTER,"@", len)) %>%
  select(`#CHROM`, start, end, name, QUAL,strand)
pbsv_vcf_lonins2 %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins.bed", col_names = FALSE)
pbsv_vcf_lonins %>%
  mutate(start = POS, end = POS, strand = "+",name = paste0(ID,"@", FILTER,"@", len)) %>% 
  filter(str_detect(sample1, "1/1")) %>% # filter out hetero insertions
  select(`#CHROM`, start, end, name, QUAL,strand) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins_homo.bed", col_names = FALSE)
pbsv_vcf_lonins %>%
  mutate(start = POS, end = POS, strand = "+",name = paste0(ID,"@", FILTER,"@", len)) %>% 
  filter(str_detect(sample1, "0/1")) %>%  # filter hetero insertions
  select(`#CHROM`, start, end, name, QUAL,strand) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins_hetero.bed", col_names = FALSE)

pbsv_vcf_londel <- pbsv_vcf %>% filter(str_detect(ID, "DEL")) %>%
  mutate(len = str_length(REF)) %>%
  filter(len >50)
pbsv_vcf_londel2 <- pbsv_vcf_londel %>%
  mutate(start = POS,
         end = POS + len-1)
pbsv_vcf_londel2 %>% mutate(name = paste0(ID,"@", FILTER ),strand = "+") %>%
  filter(str_detect(sample1, "1/1")) %>% # filter out hetero deletions
  select(`#CHROM`, start, end, name, QUAL,strand) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_del_homo.bed", col_names = FALSE)

pbsv_vcf_londel2 %>% mutate(name = paste0(ID,"@", FILTER ),strand = "+") %>%
  filter(str_detect(sample1, "0/1")) %>% # filter out hetero deletions
  select(`#CHROM`, start, end, name, QUAL,strand) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_del_hetero.bed", col_names = FALSE)

pbsv_vcf_londel2 %>% mutate(name = paste0(ID,"@", FILTER ),strand = "+") %>%
  select(`#CHROM`, start, end, name, QUAL,strand) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_del_all.bed", col_names = FALSE)

pbsv_vcf_lonins3 <- pbsv_vcf_lonins %>% rename(chr = `#CHROM`) %>%
  mutate(start = POS, end = POS+1, strand = "+",name = paste0(ID,"@", FILTER,"@", len)) %>%
  select(chr, start, end, name, QUAL,strand, ALT) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins.tsv")
pbsv_vcf_lonins3 %>%
  select(name, ALT) %>% write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins_read.tsv", col_names = FALSE) %>%
  rename(seq = ALT) %>%
  mutate(index = row_number(),name = paste0(">", name)) %>%
  gather(seq,name,key = "type",value = "fasta") %>%
  arrange(index, type) %>% select(fasta) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins_read.fasta",col_names = FALSE)

anno_pbsv <- pbsv_vcf_lonins %>%
  select(ID, FILTER,INFO,FORMAT,sample1) %>%
  separate(sample1, into = c("GT","AD","DP","SAC"),sep = ":") %>%
  mutate(certainty = if_else(str_detect(INFO,"IMPRECISE"), "imprecise","precise"),
         genotype = case_when(GT == "1/1" ~ "homo",
                              GT == "0/1" ~ "hetero",
                              TRUE ~ "notsure"),
         AD = AD %>% str_remove("[:digit:]+,") %>% as.integer(),
         DP = DP %>% as.integer()) %>%
  select(ID,FILTER,AD,DP,certainty,genotype)
anno_pbsv %>% write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins_annotation.tsv")
