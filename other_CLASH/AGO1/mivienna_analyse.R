#######################################################################################################
# making table
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)
library(gplots)
library(RColorBrewer)
library(gtools)
library(gridExtra)
library(stringi)


ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
filedir <- "/media/hermione/AGO1_CLASH/vienna/"
outdir <- "/media/hermione/AGO1_CLASH/"

############################################################
#set factor
#this is used for select read-target.
targetRNAlevel <- c("rRNA", "tRNA", "snRNA", "snoRNA", "mtrRNA", "miscRNA", "Ig", "prtr","mRNA",
                    "pseudo", "lincRNA",  "microRNA")

files3 <- dir(filedir, pattern = "_vienna.tsv$", full.names = TRUE)
files4 <- files3 %>%
  str_replace_all("_vienna.tsv$", "") %>%
  str_replace(paste0(filedir, "/"), "")
df2 <- vector("list", length(files3))
for (i in seq_along(files3)) {
  df2[[i]] <- read_tsv(files3[[i]], col_types = cols(dG = col_double()))
  df2[[i]][["sample"]] <- files4[[i]]
}
Vienna1 <- df2 %>% bind_rows() %>%
  mutate(dG = if_else(!(str_detect(miRNAmatch, "\\(")), Inf, dG))

Vienna1 %>% write_csv(paste0(outdir, "vienna1fig/clash.csv"))

Vienna_select <- Vienna1 %>%
  group_by(sample, readname) %>%
  nest_legacy() %>%
  mutate(data2 = map(data, function(x) {x %>% filter(dG == min(dG)) %>%
      filter(abs(milen-21.9)==min(abs(milen-21.9)))})) %>%
  select(-data) %>% unnest_legacy() %>%
  write_csv(paste0(outdir, "vienna1fig/clash_selection.csv"))

Vienna_sle <- Vienna1 %>%
  select(sample, readname,dG, miRNA,milen, miRNAseq, targettype,
         targetlen,targetseq,miRNAmatch, targetmatch) %>%
  group_by(sample, readname) %>%
  nest_legacy() %>%
  mutate(data2 = map(data, function(x) {x %>% filter(dG == min(dG)) %>%
      filter(abs(milen-21.9)==min(abs(milen-21.9))) %>% distinct() %>%
      mutate(targettype = targettype %>% factor(levels = targetRNAlevel)) %>%
      arrange(targettype) %>%
      slice_max(n=1, desc(targettype), with_ties = FALSE)})) %>%
  select(-data) %>% unnest_legacy()
Vienna_sle %>% write_csv(paste0(outdir, "vienna1fig/clash_colapse.csv"))
