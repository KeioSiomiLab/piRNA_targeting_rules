#
# /usr/bin/Rscript --silent --slave --vanilla clash_make_vienna_shuf2.R shuffle/${X}_chimera.tsv shuffle/ shuffle/${X}_vienna_shuffle.tsv;
#######################################################################################################
# clash_make_vienna for shuffle
#######################################################################################################
# hybdir <- "/media/pericles/CLASH/shuffle/Vienna_sle_chimera.tsv"
# rnaplexname <- "/media/pericles/CLASH/shuffle/"
# outname <- "/media/pericles/CLASH/shuffle/Vienna_sle_vienna_shuffle.tsv"
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

library(tidyverse)


hybdir <- commandArgs(trailingOnly = TRUE)[1]
rnaplexname <- commandArgs(trailingOnly = TRUE)[2]
outname <- commandArgs(trailingOnly = TRUE)[3]

hybtest <- read_tsv(hybdir)

files3 <- dir(rnaplexname, pattern = ".rnaplex$", full.names = TRUE)

rnaplex <- map(files3, function(x) {read_tsv(x, 
                  col_names = "data")}) %>% 
  bind_rows() %>% 
  mutate(index = case_when(str_detect(data,"piRNA") ~ "piRNA",
                           str_detect(data,">\\d+_\\d+") ~ "readname",
                           TRUE ~ "match"),
         counter = if_else(index =="piRNA", 1,0) %>% cumsum()) %>% 
  pivot_wider(names_from = "index",values_from = data) %>% 
  drop_na()


rnaplex2 <- rnaplex %>% mutate(match = match %>% str_replace(" \\(", "  \\(") %>% str_replace_all("   ", "  ") %>% str_replace_all("\\) 1", "\\)  1") %>% 
                                 str_replace_all("\\. 1", "\\.  1") %>% str_replace_all(" :  ", "  :  ") %>% str_remove(":  ")) %>% 
  separate(match, sep = "\\s\\s", into = c("bind", "p3_inno", "p5_inno", "dG_raw")) %>% 
  separate(bind, sep = "&", into = c("p3bind", "p5bind")) %>% 
  separate(p3_inno, sep = ",", into = c("p3_start","p3_end")) %>% 
  separate(p5_inno, sep = ",", into = c("p5_start","p5_end")) %>% 
  mutate(p3_start = p3_start %>% str_remove_all(" ") %>% as.integer(),
         p3_end = p3_end %>% str_remove_all(" ") %>% as.integer(),
         p5_start = p5_start %>% str_remove_all(" ") %>% as.integer(),
         p5_end = p5_end %>% str_remove_all(" ") %>% as.integer()) %>% 
  mutate(piRNA = piRNA %>% str_remove(">piRNA") %>% str_remove("&\\d+"), readname = readname %>%  str_remove(">"),
         dG = dG_raw %>% str_extract("\\(.+\\)") %>%  str_remove_all("\\(|\\)") %>% as.double()) %>% 
  rename(sample =piRNA)


rnaplex3 <- hybtest %>% left_join(rnaplex2,by = c("sample","readname")) %>% 
  mutate(targetlen = targetseq %>% str_length() %>% as.integer(),
         pilen = piRNAseq %>% str_length() %>% as.integer(),
         piRNAmatch = paste0(strrep(".", p3_start-1), p3bind, strrep(".", pilen - p3_end)) %>% str_replace_all("\\)", "\\("),
         targetmatch = paste0(strrep(".", p5_start-1), p5bind, strrep(".", targetlen - p5_end)) %>% str_replace_all("\\(", "\\)")) %>% 
  select(-contains("p3"), -contains("p5"), -contains("dis"), -dG_raw,-counter) %>% 
  mutate(piRNAmatch = if_else(!(str_detect(piRNAmatch, "\\(")), strrep(".", pilen),piRNAmatch),
         targetmatch = if_else(!(str_detect(piRNAmatch, "\\(")), strrep(".", targetlen),targetmatch))
rnaplex3 %>% write_tsv(outname)








