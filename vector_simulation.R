
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
ppi <- 600
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

outdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/fig6/"

#######################################################################################################
# mutation_simulation
#######################################################################################################

library(tidyverse)

#data1000 <- tibble(start = 1:971, end = 30:1000)
#stablepoint <- 1:1000
#mutpoint <- sample(stablepoint,size = 1)
#new_stable <- stablepoint[stablepoint != mutpoint]
#new_data <- data1000 %>% filter(!(start <= mutpoint & end >= mutpoint))
#row_num <- length(new_data$start)
#mutpoint <- sample(stablepoint,size = 50)
#new_stable <- stablepoint[!(stablepoint %in% mutpoint)]
#new_data <- data1000
#for (i in seq_along(mutpoint)) {
#  new_data <- new_data %>% filter(!(start <= mutpoint[[i]] & end >= mutpoint[[i]]))
#}
#row_num <- length(new_data$start)



pilen <- 26
vector_len_list <- c(500, 1000,2000,3000,4000,5000)

df4 <- df3 <- vector("list", length = length(vector_len_list))
for (vector_len_iter in seq_along(vector_len_list)){
  vector_len <- vector_len_list[[vector_len_iter]]
  stablepoint <- 1:vector_len
  data_all <- tibble(start = 1:(vector_len-pilen+1), end = pilen:vector_len)
  iter_len <- 50
  point_num <- c(1,seq(1:11)*vector_len %/% 100)
  df3 <- vector("list", length = length(point_num))
  for (point in seq_along(point_num)) {
    df2 <- vector("list", iter_len)
    for (iter in 1:iter_len) {
      mutpoint <- sample(stablepoint,size = point_num[[point]])
      new_data <- data_all
      for (i in seq_along(mutpoint)) {
        new_data <- new_data %>% filter(!(start <= mutpoint[[i]] & end > mutpoint[[i]]))
      }
      df2[[iter]] <- length(new_data$start)
    }
    df3[[point]] <- unlist(df2) %>% 
      enframe(name ="iter_num",value = "bind_num") %>%
      mutate(mut_number = point_num[[point]])
  }
  df4[[vector_len_iter]] <- bind_rows(df3) %>% 
    mutate(vector_len = vector_len )
}
df4_res <- bind_rows(df4) %>% 
  mutate(samplename = paste0(vector_len," bp") %>% factor(levels = paste0(vector_len_list," bp"))) 

df4_res_error <- df4_res %>% 
  group_by(samplename, mut_number) %>% 
  summarise(bind_mean = mean(bind_num), bind_sd = sd(bind_num)) %>% ungroup()
solve_threshold <- df4_res %>% mutate(sample = samplename) %>% 
  group_by(sample) %>% nest_legacy() %>%
  mutate(oberon = map(data, function(x) {stats::loess(x$mut_number ~ x$bind_num, span = 0.75)}),
         titania = map2(data, oberon, function(x,y) {stats::predict(y, 300)})) %>% 
  dplyr::select(-data, -oberon) %>% unnest_legacy() %>% mutate(yposi = 700  *(row_number() %>% as.integer()) %% 2)
df4_res %>% 
  ggplot(aes(x = mut_number))+
  geom_hline(yintercept=300,alpha = 0.3, color = "gray10", linetype = "dashed")+
 # geom_smooth(aes(y = bind_num, group = vector_len,color = samplename),method = "loess", alpha = 0.2,
 #              linewidth = 1, se = FALSE, fullrange = TRUE)+
  stat_smooth(geom = "line",aes(y = bind_num, group = vector_len,color = samplename),method = "loess", alpha = 0.5,
                            linewidth = 1, se = FALSE, fullrange = TRUE)+
  geom_point(aes(y = bind_num, color = samplename),alpha = 0.05)+
  geom_errorbar(data = df4_res_error,
                aes(ymin= bind_mean - bind_sd, ymax = bind_mean + bind_sd, group = vector_len),
                alpha = 0.15, color = "black",width = 10) +
  geom_segment(data = solve_threshold, aes(x=titania , xend = titania + 8,y = 300, yend =yposi * 5/7),
               color = "gray30", alpha = 0.8)+
  geom_text(data = solve_threshold, aes(x=titania + 10, label = round(titania, digits = 1),y = 50 + yposi-200), size = 4)+
#  annotate("text", x = -20, y = 300, color = "black", label = "300", size = 2.5)+
  guides(group ="none")+
  labs(x = "mutation number (bp)", y = "piRNA variety (count)",color = "binding length")+
  theme_minimal() + 
  theme(legend.position = c(0.75,0.65), panel.grid.minor = element_blank())+
  ggsave(paste0(outdir,"vector_simulation.png"), width =3, height =3, dpi = 600)+
  ggsave(paste0(outdir,"vector_simulation.pdf"), width =3, height =3)

#     text     line     
# odd   550    700
# even -100      0

 

  
