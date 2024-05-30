
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

library(tidyverse)

filedir <- "/media/ariel/Isoseq/process/"

file_flnc <- dir(filedir, pattern = "_flnc.report.csv$", full.names = TRUE)
file_flnc2 <- file_flnc %>% 
  str_replace_all("_flnc.report.csv$", "") %>%
  str_replace(paste0(filedir, "/"), "") 
file_flnc3 <- file_flnc2 %>% 
  str_remove("_")
df <- vector("list", length(file_flnc))
for (i in seq_along(file_flnc)) {
  df[[i]] <- read_csv(paste0(filedir, file_flnc2[[i]],"_flnc.report.csv"))
  df[[i]][["sample"]] <- file_flnc3[[i]]
}
data1 <- bind_rows(df) %>% select(id,sample) %>% 
  rename(primer = sample)

data1 %>% 
  write_csv(paste0(filedir, "flnc_report.csv"))


