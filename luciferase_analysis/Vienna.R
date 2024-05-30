

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/CLASH/vienna/")
suppressMessages(suppressWarnings(require(tidyverse)))

#/usr/bin/Rscript --silent --slave --vanilla piRNA_make_binding.R hyb_temp/{/.}.tsv rnaplex/{/.}.rnaplex binding/{/.};


ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}


filedir <- "/media/hermione/piRNA_luc/binding/"
outdir <- "/media/hermione/piRNA_luc/figure/"

piRNA_expression <- read_tsv("/media/pericles/CLASH/piRNA/piRNA_CPM.tsv") %>%
  select(name,seq,piRNA2,piRNA3, SRR2749801,SRR2749802,SRR9158321, CPM_avg)


files5 <- dir(filedir, pattern = "_struct.tsv$", full.names = TRUE)
files6 <- files5 %>% str_replace_all("_struct.tsv$", "") %>% str_replace(paste0(filedir, "/"), "")
df3 <- vector("list", length(files5))
for (i in seq_along(files5)) {
  df3[[i]] <- read_tsv(files5[[i]])
  df3[[i]][["sample"]] <- files6[[i]]
}

Vienna_struct <- df3 %>% bind_rows() %>%
  left_join(piRNA_expression %>% rename(piRNA =name,piRNAseq = seq), by = c("piRNA","piRNAseq"))

# functional comparison
Vienna_struct %>% mutate(picount= str_count(piRNAmatch,"\\("), pilen = str_count(piRNAmatch,"\\(|\\.") ) %>% 
  mutate(pistruct = case_when(!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)")) ~ "perfect match",
         (str_detect(piRNAmatch,"^\\.\\(\\(\\(") & pilen ==(picount + 1)) &!(str_detect(targetmatch,"\\)\\.+\\)"))  ~ "functional",
                          TRUE ~ "not functional")) %>% 
  select(target,targetseq,piRNA,piRNAseq,dG, piRNAmatch,targetmatch,ct,pistruct) %>% 
  write_csv(paste0(outdir,"mdg1_construct_struct.csv"))



