# make piRNA database
# setwd("/media/pericles/CLASH/piRNA/")
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/CLASH/piRNA/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(scales)))
suppressMessages(suppressWarnings(require(ggupset)))
suppressMessages(suppressWarnings(require(UpSetR)))
options(digits=7)

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

fasta_make <- function(file, threshold){
  #file = "SRR8791716"
  #threshold =10
  #read file
  tx1 <- read_tsv(paste0(file, ".fa"), col_names = "seq") %>%
    filter(!(str_detect(seq, "length"))) %>%
    group_by(seq) %>% summarise(count = n()) %>%
    mutate(CPM = count*1000000/sum(count))
  tx1 %>% write_tsv(paste0("otherdata/piRNA_", file, "_for_graph.tsv"))
  #filter and index
  tx3 <- tx1 %>% filter(CPM >= threshold) %>% filter(count >= 3) %>%
    mutate(annotation= paste0("piRNA", row_number()),
           index = row_number())
  tx3 %>% write_tsv(paste0("piRNA_", file, ".tsv"))
  #to fasta
  fasta <- tx3 %>% mutate(annotation2 = paste0(">", annotation)) %>%
    select(-count, -CPM) %>%
    gather(seq, annotation2, key = "type", value = "fasta") %>%
    arrange(index, type) %>% select(fasta)
  fasta %>% write_tsv(paste0("piRNA_", file, ".fa"), col_names = FALSE)
}
fasta_make("SRR2749801", 1)
fasta_make("SRR2749802", 1)
fasta_make("SRR9158321", 1)
fasta_make("piRNA1", 1)
fasta_make("piRNA2", 1)
fasta_make("piRNA3", 1)

fasta_graph_raw <- function(file){
  tx3 <- read_tsv(paste0("otherdata/piRNA_", file, "_for_graph.tsv")) %>%
    mutate(type = file)
  return(tx3)
}
piRNA3 <- fasta_graph_raw("SRR2749801")
piRNA4 <- fasta_graph_raw("SRR2749802")
piRNA5 <- fasta_graph_raw("SRR9158321")
piRNA6 <- fasta_graph_raw("piRNA2")
piRNA7 <- fasta_graph_raw("piRNA3")
piRNA9 <- fasta_graph_raw("piRNA1")

piRNAseq_ano <- tibble(type = c("SRR2749801", "SRR2749802", "SRR9158321","piRNA1","piRNA2","piRNA3"),
                       type2 = c("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#4","piRNA#5","piRNA#6"))
piRNAgraph <- bind_rows(piRNA3, piRNA4, piRNA5,piRNA6,piRNA7,piRNA9) %>%
  left_join(piRNAseq_ano, by = "type")
threshold <- 1
piRNAindex <- piRNAgraph %>% group_by(type2) %>%
  summarise(threshold = threshold * sum(count/1000000)) %>% ungroup() %>%
  mutate(threshold = if_else(threshold >= 3, threshold, 3))
piRNAgraph %>% group_by(type2, count) %>%
  summarise(count2 = n()) %>% ungroup() %>%
  ggplot() +
  geom_point(aes(x = count, y = count2, color = type2), size = 1, alpha = 0.7, shape= 19) +
  geom_vline(data = piRNAindex, aes(xintercept = threshold), color = "red") +
  geom_text(data = piRNAindex, aes(x = threshold , y = 100000, label = paste0(" threshold = ",round(threshold, digits = 1))),
            hjust=-0.2, size = 3) +
  xlab("piRNA enrichment of each data") +
  ylab("count of same enrichment") +
  scale_x_log10(breaks = 10 ^ (0:5),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = 10 ^ (0:6),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  facet_wrap( ~ type2)+
  theme_minimal()+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  guides(color = "none") +
  ggsave("figures/piRNA_threshold.png" , width =7, height =5, dpi = 300)





fasta_make2 <- function(file){
  tx3 <- read_tsv(paste0("piRNA_", file, ".tsv")) %>%
    mutate(type = file)
  return(tx3)
}
piRNA3 <- fasta_make2("SRR2749801")
piRNA4 <- fasta_make2("SRR2749802")
piRNA5 <- fasta_make2("SRR9158321")
piRNA6 <- fasta_make2("piRNA2")
piRNA7 <- fasta_make2("piRNA3")
piRNA9 <- fasta_make2("piRNA1")

piRNAseq_ano <- tibble(type = c("SRR2749801", "SRR2749802", "SRR9158321","piRNA1","piRNA2","piRNA3"),
                       type2 = c("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#4","piRNA#5","piRNA#6"))
piRNAall <- bind_rows(piRNA3, piRNA4, piRNA5,piRNA6,piRNA7,piRNA9) %>%
  left_join(piRNAseq_ano, by = "type")



piRNAcheck <- piRNAall %>% select(seq, type) %>%
  group_by(seq) %>%
  summarize(count = length(type),
            type = list(type))

ppi <- 300
png("figures/upset.png", width = 11*ppi, height = 5*ppi, res = ppi)
piRNAall %>% select(seq, type2) %>%
  distinct(seq,type2, .keep_all=TRUE) %>%
  unnest_legacy() %>%
  mutate(GenreMember=1) %>%
  pivot_wider(names_from = type2, values_from = GenreMember, values_fill = list(GenreMember = 0)) %>%
  as.data.frame() %>%
  UpSetR::upset(sets = c("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#4", "piRNA#5","piRNA#6"),sets.bar.color = "#56B4E9", order.by="freq",
                queries = list(list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#4", "piRNA#5"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#4", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#4", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#3", "piRNA#4", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#2", "piRNA#3", "piRNA#4", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#4"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#5"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#4", "piRNA#5"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#4", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#3", "piRNA#4", "piRNA#5"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#3", "piRNA#4", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#3", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#4", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               
                               list(query=intersects, params=list("piRNA#2", "piRNA#3", "piRNA#4", "piRNA#5"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#2", "piRNA#3", "piRNA#4", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#2", "piRNA#3", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#2", "piRNA#4", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#3", "piRNA#4", "piRNA#5", "piRNA#6"), color="orange", active=T),
                               list(query=intersects, params=list("piRNA#1", "piRNA#2", "piRNA#3", "piRNA#4", "piRNA#5", "piRNA#6"), color="orange", active=T)),
                nintersects = NA, mb.ratio = c(0.7, 0.3), mainbar.y.label = "count of piRNA", sets.x.label = "count of each data")
dev.off()


fasta_all <- piRNAcheck %>% filter(count >= 1) %>%
  filter(!(str_detect(seq,"^AAAAAAAAAAAAAAAAAAAAAA"))) %>% 
  select(seq) %>%
  mutate(annotation= paste0("piRNA", row_number()),
         index = row_number(),
         annotation2 = paste0(">", annotation, "_piRNA_target_piRNA")) %>%
  gather(seq, annotation2, key = "type", value = "fasta") %>%
  arrange(index, type) %>% select(fasta)
fasta_all %>% write_tsv("piRNA_all.fasta", col_names = FALSE)

