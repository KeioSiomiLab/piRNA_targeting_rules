
#######################################################################################################
# chimera processing 1
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
# /usr/bin/Rscript chimera_process1.R sam/{/.}_fusionChimeric.out.junction sam/{/.}_fusionSJ.out.tab overlap/{/.}


suppressMessages(suppressWarnings(require(tidyverse)))

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
# CJ
# chimerafile <- "/media/hermione/CAGE/fusion/siEGFPrep1_1_STARChimeric.out.junction"
# SJfile  <- "/media/hermione/CAGE/fusion/siEGFPrep1_1_STARSJ.out.tab"
# outname <- "/media/hermione/CAGE/overlap/siEGFPrep1_1"


chimerafile <- commandArgs(trailingOnly = TRUE)[1]
SJfile <- commandArgs(trailingOnly = TRUE)[2]
outname <- commandArgs(trailingOnly = TRUE)[3]

chimera_junction1 <- read_tsv(chimerafile,
                              col_types = cols(brkpt_donorA = col_integer(),
                                               brkpt_acceptorB = col_integer()),
                              comment = "#")

chimera_junction2 <- chimera_junction1 %>% filter(chr_donorA != chr_acceptorB) %>% 
  select(chr_donorA, brkpt_donorA, strand_donorA, chr_acceptorB, brkpt_acceptorB, strand_acceptorB, junction_type,read_name,num_chim_aln,PEmerged_bool) %>% 
  mutate(read1_start = brkpt_donorA-1L,
         read1_end =  brkpt_donorA,
         read2_start = brkpt_acceptorB-2L,
         read2_end = brkpt_acceptorB-1L,
         count2 = 1/num_chim_aln) %>% 
group_by(chr_donorA,read1_start,read1_end,chr_acceptorB,read2_start,read2_end,strand_donorA,strand_acceptorB) %>% 
  summarise(count = sum(count2)) %>% 
  ungroup() %>% 
  mutate(index = paste0("CJ",row_number())) %>% 
  filter(read1_start >= 0 & read2_start >= 0) %>% 
  filter(!(str_detect(chr_donorA,"chr")) | !(str_detect(chr_acceptorB,"chr")))

# write output
bind_rows(chimera_junction2 %>% select(chr_donorA,read1_start,read1_end,index,count,strand_donorA)%>%
            mutate(anno = "five_prime", type = "CJ") %>% rename(chr = chr_donorA,start = read1_start,end = read1_end,strand = strand_donorA),
          chimera_junction2 %>% select(chr_acceptorB, read2_start,read2_end,index,count,strand_acceptorB)%>%
            mutate(anno = "three_prime", type = "CJ") %>% rename(chr = chr_acceptorB,start = read2_start,end = read2_end, strand = strand_acceptorB)) %>%
  arrange(chr,start,end) %>%
  write_tsv(paste0(outname,"_junction.bed"),
            col_names = FALSE)

chimera_junction2 %>%
  write_tsv(paste0(outname,"_junction.result"))
