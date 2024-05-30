.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}


filedir <- "/media/hermione/piRNA_luc/process/"

flam1 <- read_tsv(paste0(filedir, "piRNA_query1.sam"), 
                 col_names = c("target", "flag", "target2", "start", "qual", "cigar"),
                 skip = 3) %>% arrange(start)

flam2 <- read_tsv(paste0(filedir, "piRNA_query2.sam"), 
                  col_names = c("target", "flag", "target2", "start", "qual", "cigar"),
                  skip = 3) %>% arrange(start)

# make table manually!


# 1. flam rep1
# piRNA134303
# 7139-1=7138 ~ 7138+27=7165
# piRNA268036
# 79621-1=79620 ~ 79620+25=79645
# mdg1_repeat
# 12966-1=12965 ~ 13245+28=13273

# 2. flam rep2
# piRNA134303
# 7139-1=7138 ~ 7138+27=7165
# piRNA268036
# 167499-1=167498 ~ 167498+25=167523
# 313534-1=313533 ~ 313533+25=313558
# 403830-1=403829 ~ 403829+25=403854
# mdg1_repeat
# 12967-1=12966 ~ 13218+28=13246

flam_region_anno <-tribble(
  ~name, ~sample, ~start, ~end, ~regstart, ~regend,
  "piRNA134303","flam_hap1", 7108,7195,7138,7165,
  "piRNA268036","flam_hap1", 79590,79675,79620,79645,
  "mdg1","flam_hap1", 12915,13323, 12965,13273,
  "piRNA134303","flam_hap2", 7108,7195,7138,7165,
  "piRNA268036_1","flam_hap2", 167468,167553,167498,167523,
  "piRNA268036_2","flam_hap2", 313504,313588,313533,313558,
  "piRNA268036_3","flam_hap2", 403799,403884,403829,403854,
  "mdg1","flam_hap2", 12916,13296,12966,13246)

flam_region_anno %>% 
  write_tsv(paste0(filedir, "region_info.tsv"))



