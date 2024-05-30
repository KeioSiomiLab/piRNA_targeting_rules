# TEbreak result
# /usr/bin/Rscript --silent --slave --vanilla /media/pericles/TEfind/TEbreak_res_process.R;
#####################################
# TE break
#####################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)
tebreak <- read_tsv("/media/pericles/TEfind/tebreakres/test.tsv")
tebreak2 <- tebreak %>% 
  mutate(ins_site = ((`5_Prime_End`+`3_Prime_End`) /2) %>% round(digits = 0),
         end = ins_site,
         strand = case_when(Orient_5p == "+" & Orient_3p == "+" ~ "+",
                            Orient_5p == "-" & Orient_3p == "-" ~ "-",
                            Orient_5p == "None" & Orient_3p == "None" ~ "Unknown",
                            Orient_5p == "None" ~ Orient_3p,
                            Orient_3p == "None" ~ Orient_5p,
                            TRUE ~ "both"),
         score = 1) %>% 
  select(Chromosome,ins_site,end, Superfamily, score, strand) %>% 
  arrange(Chromosome, ins_site)
tebreak2 %>% write_tsv("/media/pericles/TEfind/tebreakres/TEinssite.bed",
                       col_names = FALSE)

