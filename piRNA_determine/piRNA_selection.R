# extract 1-6 column and remove inserted/deleted/unfully mapped reads.

# /usr/bin/Rscript --silent --slave --vanilla /media/pericles/CLASH/piRNA_selection.R gene/{/.}_map.sam gene/{/.}_map2.tsv;

######################
# piRNA selection
######################

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

#hybdir <- "/media/pericles/CLASH/piRNA/gene/piRNA_all_map.sam"
#outdir <- "/media/pericles/CLASH/piRNA/gene/piRNA_all_map2"

suppressMessages(suppressWarnings(require(tidyverse)))

hybdir <- commandArgs(trailingOnly = TRUE)[1]
outdir <- commandArgs(trailingOnly = TRUE)[2]

hyb1 <- read_tsv(hybdir, 
                 col_names = c("name", "flag", "target", "start", "qual", "cigar"))

hyb2 <- hyb1 %>% 
  filter(!(str_detect(cigar, "I|D|S")))

hyb2 %>% write_tsv(paste0(outdir,".tsv"))

bed <- hyb2 %>% 
  mutate(match = cigar %>% str_remove("M") %>% as.integer(),
         start = start-1,
         end =start + match,
         target = target %>% str_extract("FBgn\\d+_") %>% str_remove("_")) %>% 
  select(target,start,end,name) %>% 
  drop_na() %>% 
  write_tsv(paste0(outdir,".bed"),col_names = FALSE)

