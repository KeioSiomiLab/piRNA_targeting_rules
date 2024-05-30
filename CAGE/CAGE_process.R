#######################################################################################################
# CAGE analysis1
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
# CAGE process
library(GGally)
library(tidyverse)
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
figoutdir <- "/media/hermione/CAGE/figures/"
ctssdir <- "/media/hermione/CAGE/ctss/"
cagedir <- "/media/hermione/CAGE/genome/"
file_cage <- dir(cagedir, pattern = ".bed$", full.names = TRUE)
file_cage2 <- file_cage  %>%
  str_replace_all(".bed$", "") %>%
  str_replace(paste0(cagedir, "/"), "")
cagedf <- vector("list", length(file_cage))
for (i in seq_along(file_cage)) {
  cagedf[[i]] <- vroom::vroom(file_cage[[i]],
               col_names = c("chr1","start1","end1","chr2","start2","end2","name","qual","strand1","strand2")) %>%
    mutate(point = if_else(strand1 == "+",start1+1, end1-1)) %>%
    group_by(chr1, point,strand1) %>% count()
  cagedf[[i]][["sample"]] <- file_cage2[[i]]
}
cage_all <- cagedf %>% bind_rows()

cage_merge_rep <- cage_all %>%
  pivot_wider(names_from = sample, values_from = n, values_fill =list(n=0)) %>% ungroup() %>%
  mutate(siEGFP_rep1 = siEGFPrep1_1+siEGFPrep1_2,
         siEGFP_rep2 = siEGFPrep2_1+siEGFPrep2_2,
         siEGFP_rep3 = siEGFPrep3_1+siEGFPrep3_2,
         siPIWI_rep1 = siPIWIrep1_1+siPIWIrep1_2,
         siPIWI_rep2 = siPIWIrep2_1+siPIWIrep2_2,
         siPIWI_rep3 = siPIWIrep3_1+siPIWIrep3_2) %>%
  select(chr1,point,strand1, siEGFP_rep1,siEGFP_rep2,siEGFP_rep3,siPIWI_rep1,siPIWI_rep2,siPIWI_rep3) %>%
  pivot_longer(-c(chr1,point, strand1), names_to = "sample", values_to = "count") %>%
  filter(count != 0)
tx <- cage_merge_rep %>% arrange(sample,chr1, point) %>%
  group_by(sample) %>% nest_legacy()
tx2 <- map(1:length(tx$sample), function(x) data.frame(tx$data[[x]]))
names(tx2) <- tx$sample
data <- map(1:length(tx2), function(x) write_tsv(tx2[[x]],
                                                 paste0(ctssdir, "cage", names(tx2[x]), ".ctss"),
                                                 col_names = FALSE))

tx <- cage_merge_rep %>% mutate(start =  point-1,
                                end =  point,
                                name = paste0("CAGE", row_number()),
                                score = 0) %>%
  select(sample, chr1,start, end, name, count,strand1) %>%
  arrange(sample,chr1, start) %>%
  group_by(sample) %>% nest_legacy()
tx2 <- map(1:length(tx$sample), function(x) data.frame(tx$data[[x]]))
names(tx2) <- tx$sample
data <- map(1:length(tx2), function(x) write_tsv(tx2[[x]],
                                                 paste0(ctssdir, "cage", names(tx2[x]), ".bed"),
                                                 col_names = FALSE))
