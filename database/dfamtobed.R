#########################################
# dfam data to bed files
#########################################
# /usr/bin/Rscript --silent --slave --vanilla dfamtobed.R database/dm6_dfam.nrph.hits database/dm6_dfam_nrph;
# dfam <- "/media/pericles/TEfind/database/dm6_dfam.nrph.hits"
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))

dfam <- commandArgs(trailingOnly = TRUE)[1]
outname <- commandArgs(trailingOnly = TRUE)[2]

dfam3 <- read_tsv(dfam) 
dfam4 <- dfam3 %>%
  rename(seq_name = `#seq_name`) %>% 
  rename(chr = seq_name,
    name = family_name,
    start = `ali-st`,
    end = `ali-en`,
    score = kimura_div) %>% 
  mutate(new_start = if_else(strand =="+", start, end) %>% as.integer(),
         new_end = if_else(strand == "+", end, start) %>% as.integer()) %>% 
  select(chr, new_start, new_end, name, score, strand)

dfam4 %>% write_tsv(paste0(outname, "_full.bed"), col_names = FALSE)
dfam4 %>% filter(chr %in% paste0("chr", c("2L","2R","3L","3R","4","X","Y"))) %>% 
  write_tsv(paste0(outname, ".bed"), col_names = FALSE)

# add class
dfamclass <- read_tsv("/media/pericles/TEfind/database/class_process_all_species.tsv")
dfam5 <- dfam4 %>% left_join(dfamclass %>% rename(name = TE), by = "name") %>%
  mutate(subfamily = if_else(is.na(subfamily), "non",subfamily)) %>% 
  unite(family,subfamily,sep = "@",col = "family")

dfam5 %>% select(chr,new_start, new_end, family,name, strand) %>% 
  write_tsv("/media/pericles/TEfind/database/annotateTE/TEannotate_tebreak.txt",col_names = FALSE)
#########################################################################################
# telomere and centromere
#########################################################################################

centel <- read_tsv("/media/pericles/TEfind/database/cytoBand.txt",
                   col_names = c("chr","start","end","dis","name")) %>% 
  filter(chr %in% paste0("chr", c("2L","2R","3L","3R","4","X","Y"))) %>% 
  filter(name == "gneg") %>% 
  arrange(chr,start,end) %>% 
  select(chr,start,end, name)
centel %>% write_tsv("/media/pericles/TEfind/database/centromere_telomere.txt",
                     col_names = FALSE)

centel <- read_tsv("/media/pericles/TEfind/database/cytoBand.txt",
                   col_names = c("chr","start","end","dis","name"),
                   col_types = cols(chr = col_character(),
                                    start = col_double(),
                                    end = col_double(),
                                    dis = col_character(),
                                    name = col_character())) %>% 
  filter(chr %in% paste0("chr", c("2L","2R","3L","3R","4","X","Y"))) %>% 
  write_tsv("/media/pericles/TEfind/database/cytoBand_igv.txt",
            col_names = FALSE)
  


