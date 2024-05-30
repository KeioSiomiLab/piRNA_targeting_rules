# make vienna table
# /usr/bin/Rscript --silent --slave --vanilla vcf_process.R;
#######################################################################################################
# vienna table make
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/hermione/cCLIP/vienna/")
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
filedir <- "/media/hermione/cCLIP/vienna/"
splicedir <- "/media/hermione/cCLIP/splice/"
outdir <- "/media/hermione/cCLIP/"

anotation1 <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_use.tsv") %>%
  select(type, name, genesymbol) %>%
  filter(type %in% c("rRNA", "tRNA", "snRNA", "snoRNA", "mRNA", "pseudogene", "ncRNA", "pre_miRNA", "miRNA")) %>%
  distinct() %>% rename(target = name,	symbol = genesymbol) %>% select(target, symbol, type)

gene_all <- read_tsv("/media/pericles/CLASH/database/dm6_gene.tsv",
                     col_names = c("target", "seq")) %>%
  separate(target, into = c("target","dis6","dis7"), sep = "_") %>% select(-dis6,-dis7) %>%
  mutate(target = target %>% str_remove("\\-"))
#######################################################################################################
# detect overlap
#######################################################################################################

typeanno <- "/media/hermione/cCLIP/typeanno/"
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
#######################################################################################################
# load vienna
#######################################################################################################
#set factor
#this is used for select read-target.
targetRNAlevel <- c("rRNA", "tRNA", "snRNA", "snoRNA", "ensemblTE", "mRNA(TE)","mRNA",
                    "pseudogene", "ncRNA", "pre_miRNA", "miRNA", NA)

files3 <- dir(filedir, pattern = "_vienna.tsv$", full.names = TRUE)
files4 <- files3 %>%
  str_replace_all("_vienna.tsv$", "") %>%
  str_replace(paste0(filedir, "/"), "")
df2 <- vector("list", length(files3))
for (i in seq_along(files3)) {
  df2[[i]] <- read_tsv(files3[[i]], col_types = cols(dG = col_double())) %>% mutate(des = "bowtie2")
  df2[[i]][["sample"]] <- files4[[i]]
}
Vienna1 <- df2 %>% bind_rows() %>%
  mutate(dG = if_else(!(str_detect(piRNAmatch, "\\(")), Inf, dG)) %>%
  left_join(anotation1, by = "target") %>%
  left_join(typeall2, by = c("sample","readname","target","bedstart","bedend")) %>%
  mutate(symbol = if_else(targettype=="ensemblTE", target, symbol),
         targettype = if_else(targettype=="ensemblTE", targettype,
                              if_else(targettype=="pbsvTE" | targettype=="rmblastTE", targettype, type)),
         targettype = if_else(targettype=="mRNA" & str_detect(targettype2,"TE_sense|TE_anti"), "mRNA(TE)", targettype)) %>%
  select(-type) %>% left_join(gene_all,by = "target") %>%
  mutate(seq = str_sub(seq,start = bedstart,end=bedstart) %>% str_to_upper()) %>%
  mutate(group = "ovary")

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
  select(sample, readname,dG, piRNA,pilen, piRNAseq, targettype,
         targetlen,targetseq,piRNAmatch, targetmatch,des,direction) %>%
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

Vienna_select <- Vienna1 %>%
  group_by(sample, readname) %>%
  nest_legacy() %>%
  mutate(data2 = map(data, function(x) {x %>% filter(dG == min(dG)) %>%
      filter(abs(pilen-26.1)==min(abs(pilen-26.1)))})) %>%
  select(-data) %>% unnest_legacy() %>%
  write_csv(paste0(outdir, "vienna1fig/clash_selection.csv"))

Vienna_sle <- Vienna1 %>%
  select(sample, readname,dG, piRNA,pilen, piRNAseq, targettype,
         targetlen,targetseq,piRNAmatch, targetmatch,des,direction) %>%
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
