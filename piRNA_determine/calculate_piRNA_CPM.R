# calculate CPM of piRNA
# This calculation includes some mistake. this mistake is corrected later.
##########################################################
#CPM cal
##########################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(GGally)))
suppressMessages(suppressWarnings(require(ggthemes)))
suppressMessages(suppressWarnings(require(ggseqlogo)))

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

filedir <- "/media/pericles/CLASH/piRNA/"
outdir <- "/media/pericles/CLASH/piRNA/"
figoutdir <- "/media/pericles/CLASH/piRNA/figout/"
piRNA <- read_tsv("/media/pericles/CLASH/piRNA/piRNA_all.tsv",
                  col_names = c("name", "seq")) %>%
  mutate(name = name %>% str_extract("piRNA\\d+"),
         pilen = str_length(seq))


files3 <- dir(filedir, pattern = ".tsv$", full.names = TRUE) %>%
  str_subset("piRNA_") %>% str_subset("piRNA_all|piRNA_CPM", negate = TRUE)
files4 <- files3 %>%
  str_replace_all(".tsv$", "") %>%
  str_replace(paste0(filedir, "/"), "") %>%
  str_remove("piRNA_")
df2 <- vector("list", length(files3))
for (i in seq_along(files3)) {
  df2[[i]] <- piRNA %>% left_join(read_tsv(files3[[i]]) %>% select(seq, count,CPM), by = "seq")
  df2[[i]][["sample"]] <- files4[[i]]
}
tx <- df2 %>% bind_rows() %>%
  mutate(CPM = if_else(is.na(CPM), 0, CPM)) %>%
  select(name, seq, CPM, sample, pilen) %>%
  pivot_wider(names_from = sample, values_from = CPM)
tx2 <- tx %>% mutate(CPM_avg = (SRR2749801 + SRR2749802 + SRR9158321 + piRNA1 + piRNA2 + piRNA3) / 6)
tx3 <- select(tx2,-seq,-name,-pilen,-CPM_avg) %>%
  rename("piRNA#1" = "SRR2749801",
         "piRNA#2" = "SRR2749802",
         "piRNA#3" = "SRR9158321",
         "piRNA#4" = "piRNA1",
         "piRNA#5" = "piRNA2",
         "piRNA#6" = "piRNA3") %>%
  mutate_at(vars(contains("piRNA")),list(~log10(.+1)))

ppi <- 300
cor_func <- function(data,mapping,method,symbol){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x,y, method = method,use = 'complete.obs')
  colFn <- colorRampPalette(c("brown1","white","dodgerblue"),interpolate = "spline")
  fill <- colFn(100)[findInterval(corr,seq(-1,1,length = 100))]
  ggally_text(label = paste(symbol,as.character(round(corr, 2))),
              mapping = aes(),
              xP = 0.5,yP = 0.5,color = "black")+
    theme(panel.background = element_rect(fill = fill))
}
theme_set(theme_minimal())
png(paste0(figoutdir,"ggally.png"), width = 5*ppi, height = 5*ppi, res = ppi)
ggpairs(data = tx3 %>% select("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#4", "piRNA#5","piRNA#6"), xlab = "log10(CPM + 1)",ylab = "log10(CPM + 1)",
        lower = list(continuous = wrap("points", alpha = 0.2, size = 0.1)),
        upper = list(continuous = wrap(cor_func, method = "spearman", symbol = expression('\u03C1 ='))))
dev.off()

tx2 %>%
  write_tsv(paste0(outdir,"piRNA_CPM.tsv"))

