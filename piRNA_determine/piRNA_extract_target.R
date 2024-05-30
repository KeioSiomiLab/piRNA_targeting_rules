# add annotation to piRNA
# /usr/bin/Rscript --silent --slave --vanilla /media/pericles/CLASH/piRNA_extract_target.R gene/{/.}_map2.tsv \
# /media/pericles/TEfind/tldrres/allTE.tsv \
# {/.}.tsv \
# anno/{/.} \
# figures/{/.};
######################
# gene_TE
######################
# sam <- "/media/pericles/CLASH/piRNA/gene/piRNA_all_map2"
# TE <- "/media/pericles/TEfind/tldrres/allTE.tsv"
# piRNAname <- "/media/pericles/CLASH/piRNA/piRNA_all.tsv"
# out <- "/media/pericles/CLASH/piRNA/anno/piRNA_all"
# figout <- "/media/pericles/CLASH/piRNA/figures/piRNA_all"
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")


suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(stringi)))
suppressMessages(suppressWarnings(require(ggseqlogo)))

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

sam <- commandArgs(trailingOnly = TRUE)[1]
TE <- commandArgs(trailingOnly = TRUE)[2]
piRNAname <- commandArgs(trailingOnly = TRUE)[3]
out <- commandArgs(trailingOnly = TRUE)[4]
figout <- commandArgs(trailingOnly = TRUE)[5]

samgene <- read_tsv(paste0(sam,".tsv"))
samgene2 <- samgene %>%
  mutate(name = name %>% str_extract("piRNA\\d+")) %>%
  separate(target, sep= "_", into = c("FBgn", "target","distype"))



overlap <- read_tsv(paste0(sam %>% str_remove("_map2"),"_overlap.tsv"),
                    col_names = c("target","bedstart","bedend","readname","gene","start","end","name","anno"))
overlap_exist <- read_tsv(paste0(sam %>% str_remove("_map2"),"_overlap2.tsv"),
                          col_names = c("target","bedstart","bedend","readname","gene","start","end","name","anno")) %>%
  mutate(exist = "yes")


overlap2 <- overlap %>% left_join(overlap_exist, by = c("target","bedstart","bedend","readname","gene","start","end","name","anno")) %>%
  mutate(exist = if_else(str_detect(anno, "TE") & is.na(exist), "no",exist),
         anno = if_else(is.na(exist),anno,paste0(anno,"_",exist))) %>%
  filter(!(anno %in% c("TE_anti_no", "TE_sense_no"))) %>% 
  select( target,bedstart,bedend,readname, name,anno) %>%
  group_by(readname,target,bedstart,bedend) %>% nest_legacy() %>%
  mutate(annoname = map(data, function(x) {str_c(x$name,collapse = "|")}),
         targettype2 = map(data, function(x) {x$anno %>% unique() %>% str_c(collapse = "|")})) %>%
  select(-data) %>%
  unnest_legacy() %>% ungroup()

samgene2_ove <- samgene2 %>% mutate(bedstart = start -1) %>%
  left_join(overlap2 %>% mutate(name = readname %>% str_extract("piRNA\\d+")) %>% rename(FBgn = target),
            by = c("name","FBgn","bedstart")) %>%
  mutate(targettype2 = if_else(is.na(targettype2),"TE",targettype2)) %>%
  filter(!(str_detect(targettype2, "TE_anti_no|TE_sense_no")))


samgene2_ove %>%
  write_tsv(paste0(out,"_annotation_raw.tsv"))

samgene3 <- samgene2_ove %>% select(name, target) %>% distinct() %>%
  group_by(name) %>% nest_legacy() %>%
  mutate(target = map(data, function(x) {str_c(x$target, collapse = "@")})) %>%
  select(-data) %>%
  unnest_legacy()
#TE name list
TEname <- read_tsv(TE) %>%
  select(name2) %>%
  mutate(TEname2 = name2 %>% str_replace("-", "|") %>% str_replace("_", "|"))

samgene4 <- samgene3 %>%
  mutate(type = case_when(str_detect(target, "rRNA|mtRNA|snRNA|snoRNA|tRNA") ~ "rmstpiRNA",
                          str_detect(target, "lncRNA:flam") ~"flampiRNA",
                          str_detect(target, "lncRNA:20A") ~"20ApiRNA",
                          str_detect(target, TEname$TEname2) ~ "TEpiRNA",
                          str_detect(target, "pbsv|rmblast") ~ "TEpiRNA2",
                          TRUE ~ "genicpiRNA"),
         titania = "map")

piRNA <- read_tsv(piRNAname,
                  col_names = c("name", "seq")) %>%
  mutate(name = name %>% str_extract("piRNA\\d+_") %>% str_remove("_")) %>%
  left_join(samgene4, by = "name") %>%
  mutate(titania = replace_na(titania, "non"),
         target = if_else(titania == "non", "nonorigin", target),
         type = if_else(titania != "map", "piRNA", type),
         name3 = paste0(name, "_", "target", "_",type))
piRNA %>% select(name, target, type) %>%
  write_tsv(paste0(out,"_annotation.tsv"))

piRNA2 <- piRNA %>% select(name3, seq) %>%
  mutate(index = row_number(),
         annotation2 = paste0(">", name3)) %>%
  gather(seq, annotation2, key = "type", value = "fasta") %>%
  arrange(index, type) %>% select(fasta)
piRNA2 %>% write_tsv(paste0(out, "2.fasta"), col_names = FALSE)

# gene region for last bases
gettype <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_RNAtype.tsv")
piRNA <- read_tsv(piRNAname, col_names = c("name", "seq"))

samgene4 <- piRNA %>% left_join(
  samgene %>% mutate(gene_name = target %>% str_extract("FBgn\\d+_") %>% str_remove("_")) %>%
    left_join(gettype %>% rename(gene_name = name), by = "gene_name") %>%
    mutate(type = if_else(is.na(type),"TE",type)),
  by ="name") %>%
  mutate(type = if_else(is.na(type),"Not-mapped",type))


gene_len <- read_tsv("/media/pericles/CLASH/database/dm6_gene.fasta.fai",
                     col_names = c("target","gene_len", "dis1","dis2","dis3")) %>% select(-contains("dis")) 
region_mapped <- samgene4 %>% select(target,start, cigar, name,seq,gene_name,type, genesymbol) %>% drop_na() %>% 
  mutate(length = str_length(seq) %>% as.integer(),
         mapped_len = cigar %>% str_remove("M") %>% as.integer()) %>% 
  left_join(gene_len, by = "target") %>% 
  filter(length == mapped_len) %>% 
  mutate(last_start = start + length -5,
         last_end = last_start + 9) %>% 
  filter(last_end <= gene_len)


figout2 <- "/media/pericles/CLASH/piRNA/gene_reg/"
out2 <- out %>% str_remove("anno/")
region_mapped %>% 
  write_tsv(paste0(figout2,out2,"_region_map.tsv"))

region_mapped %>% 
  select(target, last_start, last_end, name) %>% 
  write_tsv(paste0(figout2,out2,"_gene.bed"), col_names = FALSE)
