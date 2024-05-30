# make vienna table
# /usr/bin/Rscript  vcf_process.R;
#######################################################################################################
# vienna table make
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/CLASH/vienna/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gplots)))
suppressMessages(suppressWarnings(require(RColorBrewer)))
suppressMessages(suppressWarnings(require(gtools)))
suppressMessages(suppressWarnings(require(gridExtra)))
suppressMessages(suppressWarnings(require(stringi)))
suppressMessages(suppressWarnings(require(ggseqlogo)))
suppressMessages(suppressWarnings(require(rlang)))

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}


filedir <- "/media/pericles/CLASH/vienna/"
splicedir <- "/media/pericles/CLASH/splice/"
outdir <- "/media/pericles/CLASH/"


flamanno <-tribble(
  ~target, ~symbol, ~type,
  "FBgn9999996", "lncRNA:flamhap1","ncRNA",
  "FBgn9999998", "lncRNA:20Ahap1", "ncRNA",
  "FBgn9999997", "lncRNA:flamhap2", "ncRNA",
  "FBgn9999999", "lncRNA:20Ahap2","ncRNA")

anotation1 <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_use.tsv") %>%
  select(type, name, genesymbol) %>%
  filter(type %in% c("rRNA", "tRNA", "snRNA", "snoRNA", "mRNA", "pseudogene", "ncRNA", "pre_miRNA", "miRNA")) %>%
  distinct() %>% rename(target = name,	symbol = genesymbol) %>% select(target, symbol, type) %>% 
  bind_rows(flamanno)


piRNAfrom <- read_tsv("/media/pericles/CLASH/piRNA/anno/piRNA_all_annotation.tsv") %>%
  rename(piRNA = name,
         piRNAfrom = target) %>% select(-type)
# piRNA_ext
ori_filter <-  read_tsv("/media/pericles/CLASH/piRNA/anno/piRNA_all_annotation_raw.tsv") %>% 
  select(name,FBgn, target) %>% distinct() %>% 
  mutate(type = case_when(str_detect(target, "rRNA|mtRNA|snRNA|snoRNA|tRNA") ~ "rmstpiRNA",
                          TRUE ~ "others")) %>% 
  filter(type != "others") %>% 
  select(name,FBgn) %>% rename(piRNA= name, target = FBgn)



#######################################################################################################
# detect overlap
#######################################################################################################

typeanno <- "/media/pericles/CLASH/typeanno/"
file_typeanno <- dir(typeanno, pattern = "_overlap.tsv$", full.names = TRUE)
file_typeanno2 <- file_typeanno %>%
  str_replace_all("_overlap.tsv$", "") %>%
  str_replace(paste0(typeanno, "/"), "")
typeannodf <- vector("list", length(file_typeanno))
for (i in seq_along(file_typeanno)) {
  typeannodf[[i]] <- read_tsv(file_typeanno[[i]],
                              col_names = c("target","bedstart","bedend","readname","symbol","type",
                                            "gene","start","end","name","anno"))
  typeannodf[[i]][["sample"]] <- file_typeanno2[[i]]
}
typeall <- typeannodf %>% bind_rows()
typeall2 <- typeall %>% select(sample, target,bedstart,bedend,readname, name,anno) %>%
  group_by(sample,readname,target,bedstart,bedend) %>% nest_legacy() %>%
  mutate(annoname = map(data, function(x) {str_c(x$name,collapse = "|")}),
         targettype2 = map(data, function(x) {x$anno %>% unique() %>% str_c(collapse = "|")})) %>%
  select(-data) %>%
  unnest_legacy() %>% ungroup()
typeall2 %>% write_tsv(paste0(typeanno, "type_all.tsv"))

# pbsv
file_typeanno3 <- dir(typeanno, pattern = "_overlap2.tsv$", full.names = TRUE)
file_typeanno4 <- file_typeanno3 %>%
  str_replace_all("_overlap2.tsv$", "") %>%
  str_replace(paste0(typeanno, "/"), "")
typeannodf2 <- vector("list", length(file_typeanno3))
for (i in seq_along(file_typeanno3)) {
  typeannodf2[[i]] <- read_tsv(file_typeanno3[[i]],
                               col_names = c("target","bedstart","bedend","readname","symbol","type",
                                             "gene","start","end","name","score","strand"))
  typeannodf2[[i]][["sample"]] <- file_typeanno4[[i]]
}
typeall3 <- typeannodf2 %>% bind_rows()

typeall4 <- typeall3 %>% select(sample, readname,target,bedstart,bedend, name,strand) %>% distinct() %>%
  mutate(strand = if_else(strand == "+","TEplus","TEminus")) %>%
  group_by(sample,readname,target,bedstart,bedend) %>%
  mutate(annoname = str_c(name,collapse = "|"),
         targetype2 = str_c(strand,collapse = "|")) %>%
  ungroup() %>% select(-name, -strand) %>% distinct()
typeall4  %>%  write_tsv(paste0(typeanno, "type_all_pbsv.tsv"))

#######################################################################################################
# load vienna
#######################################################################################################
#set factor
#this is used for select read-target.
targetRNAlevel <- c("rRNA", "tRNA", "snRNA", "snoRNA", "ensemblTE","pbsvTE","rmblastTE", "mRNA(TE)","mRNA",
                    "pseudogene", "ncRNA", "pre_miRNA", "miRNA", NA)

files3 <- dir(filedir, pattern = "_vienna2.tsv$", full.names = TRUE)
files4 <- files3 %>% str_replace_all("_vienna2.tsv$", "") %>% str_replace(paste0(filedir, "/"), "")
sp1 <- dir(splicedir, pattern = "_vienna2.tsv$", full.names = TRUE)
sp2 <- sp1 %>%
  str_replace_all("_vienna2.tsv$", "") %>%
  str_replace(paste0(splicedir, "/"), "")
df2 <- vector("list", length(files3))
for (i in seq_along(files3)) {
  df2[[2*i-1]] <- read_tsv(files3[[i]], col_types = cols(dG = col_double())) %>% mutate(des = "bowtie2")
  df2[[2*i-1]][["sample"]] <- files4[[i]]
  df2[[2*i]] <- read_tsv(sp1[[i]], col_types = cols(dG = col_double())) %>% mutate(des = "hisat2")
  df2[[2*i]][["sample"]] <- sp2[[i]]
}


Vienna1_raw <- df2 %>% bind_rows() %>% filter(targetname != "mt:ori") %>% filter(sample != "rep6") %>% 
  mutate(group = if_else(str_detect(sample, "HP1"),"siHP1a", "WT"),
         dG = if_else(!(str_detect(piRNAmatch, "\\(")), Inf, dG)) %>%
  left_join(anotation1, by = "target") %>%
  left_join(piRNAfrom, by = "piRNA") %>%
  left_join(bind_rows(typeall2,typeall4), by = c("sample","readname","target","bedstart","bedend")) %>% 
  mutate(targettype2 = if_else(des == "hisat2", "CDS",targettype2)) %>% 
  mutate(symbol = if_else(targettype=="ensemblTE", target, symbol),
         targettype = case_when(targettype=="pbsvTE" | targettype=="rmblastTE" | targettype=="ensemblTE" ~ targettype,
                                TRUE ~ type)) %>% 
  mutate(targettype = if_else(targettype=="mRNA" & str_detect(targettype2,"TE_sense|TE_anti"), "mRNA(TE)", targettype)) %>% 
  select(-type) %>% 
  mutate(targettype = if_else(targetname == "mod(mdg4)", "mRNA",targettype)) %>% 
  filter(str_count(targetseq, "A")/str_length(targetseq) <= 0.9) %>% 
  filter(str_count(targetseq, "T")/str_length(targetseq) <= 0.9) %>% 
  filter(str_count(targetseq, "G")/str_length(targetseq) <= 0.9) %>% 
  filter(str_count(targetseq, "C")/str_length(targetseq) <= 0.9) %>% 
  filter(str_count(piRNAseq, "A")/str_length(piRNAseq) <= 0.9) %>% 
  filter(str_count(piRNAseq, "T")/str_length(piRNAseq) <= 0.9) %>% 
  filter(str_count(piRNAseq, "G")/str_length(piRNAseq) <= 0.9) %>% 
  filter(str_count(piRNAseq, "C")/str_length(piRNAseq) <= 0.9) 

Vienna1_raw %>% write_csv(paste0(outdir, "vienna1fig/clash_raw.csv"))

Vienna1_njdnjg <- df2 %>% bind_rows() %>%
  filter(!(str_detect(piRNAmatch, "\\(")))


same_ori <- ori_filter %>% left_join(Vienna1_raw %>% 
  select(readname, piRNA,target,targetname, targettype,direction,des,index,dG,sample, group,piRNAfrom), by = c("piRNA","target")) %>% 
  drop_na()
Vienna1 <- Vienna1_raw %>% 
  anti_join(same_ori, by = c("sample","readname","direction","des"))

Vienna1 %>% write_csv(paste0(outdir, "vienna1fig/clash.csv"))


#######################################################################################################
# stats of pi-tag/tag-pi
#######################################################################################################

Vienna_select_pitag <- Vienna1 %>%
  group_by(sample, readname) %>%
  nest_legacy() %>%
  mutate(data2 = map(data, function(x) {x %>% filter(dG == min(dG)) %>%
      filter(abs(pilen-26.1)==min(abs(pilen-26.1)))})) %>%
  select(-data) %>% unnest_legacy() %>%
  write_csv(paste0(outdir, "vienna1fig/clash_selection_all.csv"))

Vienna_sle_pitag <- Vienna1 %>%
  select(sample, readname,dG, piRNA,piRNAtype,pilen, piRNAseq, targettype,
         targetlen,targetseq,piRNAmatch, targetmatch,group,des,direction) %>%
  group_by(sample, readname) %>%
  nest_legacy() %>%
  mutate(data2 = map(data, function(x) {x %>% filter(dG == min(dG)) %>%
      filter(abs(pilen-26.1)==min(abs(pilen-26.1))) %>%
      distinct() %>%
      mutate(targettype = targettype %>% factor(levels = targetRNAlevel)) %>%
      arrange(targettype) %>%
      slice_max(n=1, desc(targettype), with_ties = FALSE)})) %>%
  select(-data) %>% unnest_legacy()
Vienna_sle_pitag %>% write_csv(paste0(outdir, "vienna1fig/clash_colapse_all.csv"))

#######################################################################################################
# save
#######################################################################################################

Vienna_select <- Vienna1 %>% filter(direction=="pi_tag") %>%
  group_by(sample, readname) %>%
  nest_legacy() %>%
  mutate(data2 = map(data, function(x) {x %>% filter(dG == min(dG)) %>%
      filter(abs(pilen-26.1)==min(abs(pilen-26.1)))})) %>%
  select(-data) %>% unnest_legacy() %>%
  write_csv(paste0(outdir, "vienna1fig/clash_selection.csv"))

Vienna_sle <- Vienna1 %>% filter(direction=="pi_tag") %>%
  select(sample, readname,dG, piRNA,piRNAtype,pilen, piRNAseq, targettype,
         targetlen,targetseq,piRNAmatch, targetmatch,group,des,direction) %>%
  group_by(sample, readname) %>%
  nest_legacy() %>%
  mutate(data2 = map(data, function(x) {x %>% filter(dG == min(dG)) %>%
      filter(abs(pilen-26.1)==min(abs(pilen-26.1))) %>%
      distinct() %>%
      mutate(targettype = targettype %>% factor(levels = targetRNAlevel)) %>%
      arrange(targettype) %>%
      slice_max(n=1, desc(targettype), with_ties = FALSE)})) %>%
  select(-data) %>% unnest_legacy()
Vienna_sle %>% write_csv(paste0(outdir, "vienna1fig/clash_colapse.csv"))

################################################
#for shuffle seq
################################################
filedir <- "/media/pericles/CLASH/vienna/"
shufdir <- "/media/pericles/CLASH/shuffle/"
outdir <- "/media/pericles/CLASH/"
Vienna_sle <- read_csv(paste0(outdir, "vienna1fig/clash_colapse.csv"), col_types = cols(dG = col_double()))

tag_shuf <- Vienna_sle %>% select(sample, readname, piRNAseq,targetseq) %>%
  mutate(targetseq = targetseq %>% stri_rand_shuffle(),
         index = row_number()) %>%
  write_tsv(paste0(shufdir, "Vienna_sle_chimera.tsv"))
tag_shuf2 <- tag_shuf %>% select(sample, index, readname, piRNAseq,targetseq) %>%
  mutate(piRNA = paste0(">piRNA",sample,"&",row_number()), readname = paste0(">",readname),index = row_number()) %>%
  pivot_longer(cols = c("piRNA","piRNAseq","readname","targetseq"), names_to = "type", values_to = "fasta") %>%
  select(sample,fasta)
tx <- tag_shuf2 %>% group_by(sample) %>% nest_legacy()
tx2 <- map(1:length(tx$sample), function(x) data.frame(tx$data[[x]]))
names(tx2) <- tx$sample
data <- map(1:length(tx2), function(x) write_tsv(tx2[[x]],
                                                 paste0(shufdir, "Vienna_sle_",names(tx2[x]),"_rnaplex.fasta"),
                                                 col_names = FALSE))
