
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)
#library(coin)
library(rstatix)
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
library(GGally)

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




##########################################
# output analysis 
##########################################

read_count1 <- read_tsv("/media/hermione/CAGE/comp_stats.tsv")
#samplename <- c("rep1", "rep2","rep3","rep4","rep5","rep6","rep7","rep8")
read_count2 <- read_count1 %>%
  filter(!(str_detect(file,"_2_trim.fastq.gz"))) %>% 
  mutate(sample = file %>% str_remove("_1_trim.fastq.gz") %>% str_remove("trim\\/")) %>% 
  mutate(tes_rep = sample %>% str_extract("_\\d") %>% str_remove("_") %>% as.integer(),
         sample = sample %>% str_remove("_\\d")) %>%
  group_by(sample) %>% 
  summarise(num_seqs = sum(num_seqs)) %>% ungroup() %>% 
  select(sample, num_seqs)

output <- "/media/hermione/CAGE/fusion_figure/"

chimeradir <- "/media/hermione/CAGE/chimera_res/"
chimera_file <- dir(chimeradir, pattern = "_chimera.tsv$", full.names = TRUE)
chimera_file2 <- chimera_file %>% str_remove("_chimera.tsv$") %>%
  str_remove(paste0(chimeradir, "/"))
chimeradf <- vector("list", length(chimera_file))
for (i in seq_along(chimera_file)) {
  chimeradf[[i]] <- read_tsv(chimera_file[[i]])
  chimeradf[[i]][["sample"]] <- chimera_file2[[i]]
}

chimeraall <- chimeradf %>% bind_rows() %>% 
  mutate(fusionclass = case_when(!str_detect(chr_donorA, "chr") & !str_detect(chr_acceptorB, "chr") ~ "TE_TE", TRUE ~ "valid"),
         TEtype = case_when(!str_detect(chr_donorA, "chr")~ chr_donorA, !str_detect(chr_acceptorB, "chr")~ chr_acceptorB, TRUE ~ "negative")) %>% 
  mutate(tes_rep = sample %>% str_extract("_\\d") %>% str_remove("_") %>% as.integer(),
         sample = sample %>% str_remove("_\\d")) %>% 
  group_by(sample, symbol, chr_donorA, read1_start, read1_end, chr_acceptorB,read2_start, read2_end,strand_donorA, strand_acceptorB,  gene_name, strand, rnatype, fusionclass,TEtype) %>%
  summarise(count = sum(count)) %>% ungroup() %>% 
  left_join(read_count2, by = "sample") %>% 
  mutate(CPM = count*1000000/num_seqs) %>% 
  select(-num_seqs)

chimeraall2 <- chimeraall %>% 
  select(-count) %>% 
  pivot_wider(names_from = "sample",values_from = "CPM", values_fill =list(CPM=0))

chimera_write <- chimeradf %>% bind_rows() %>% 
  mutate(fusionclass = case_when(!str_detect(chr_donorA, "chr") & !str_detect(chr_acceptorB, "chr") ~ "TE_TE", TRUE ~ "valid"),
         TEtype = case_when(!str_detect(chr_donorA, "chr")~ chr_donorA, !str_detect(chr_acceptorB, "chr")~ chr_acceptorB, TRUE ~ "negative")) %>% 
  mutate(tes_rep = sample %>% str_extract("_\\d") %>% str_remove("_") %>% as.integer(),
         sample = sample %>% str_remove("_\\d")) %>% 
  group_by(sample, symbol, chr_donorA, read1_start, read1_end, chr_acceptorB,read2_start, read2_end,strand_donorA, strand_acceptorB,  gene_name, strand, rnatype, fusionclass,TEtype) %>%
  summarise(count = sum(count)) %>% ungroup() %>% 
  left_join(read_count2, by = "sample") %>% 
  filter(fusionclass == "valid") %>%
  write_csv(paste0(output,"all_chimera_list.csv"))

graph <- ggpairs(data = chimeraall2 %>% filter(fusionclass == "valid") %>% select(contains("rep")) %>%  mutate(across(everything(), ~ log2(.x +1))), xlab = "log2(CPM + 1)",ylab = "log2(CPM + 1)",
                 lower = list(continuous = wrap("points", alpha = 0.2, size = 0.1)),
                 upper = list(continuous = wrap(cor_func, method = "pearson", symbol = expression('\u03C1 ='))))
png(paste0(output,"ggally_fusion_all.png"), width = 5*ppi, height = 5*ppi, res = ppi)
print(graph)
dev.off()


tes2 <- chimeraall2 %>% filter(fusionclass == "valid") %>% 
  mutate(type = case_when(str_detect(TEtype, "blood") ~ "dependent",
                          str_detect(TEtype, "mdg1") ~ "dependent",
                          str_detect(TEtype, "297") ~ "dependent",
                          str_detect(TEtype, "412") ~ "dependent",
                          str_detect(TEtype, "gypsy") ~ "dependent",
                          str_detect(TEtype, "gtwin") ~ "dependent",
                          str_detect(TEtype, "17.6") ~ "dependent",
                          str_detect(TEtype, "Quasimodo") ~ "dependent",
                          str_detect(TEtype, "Stalker") ~ "dependent",
                          TRUE ~ "independent"),
         siEGFP =(siEGFPrep1+siEGFPrep2+siEGFPrep3)/3,
         siPIWI =(siPIWIrep1+siPIWIrep2+siPIWIrep3)/3,
         logFC = log2((siPIWI+1)/(siEGFP +1)))%>% 
  select(-contains("rep"))
tes2 %>% 
  ggplot(aes(x = log2(siEGFP+1),y = log2(siPIWI+1), color = type))+
  geom_point(size = 1)+
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_wrap(~type) +
  theme_minimal()+
  guides(color ="none")+
  ggsave(paste0(output,"EGFP_PIWI_scatter2",".png"), width =4.5, height =3.5, dpi = 300)
#######################



# http://www.sthda.com/english/wiki/paired-samples-wilcoxon-test-in-r
cal_box_TE <- function(data, TE ){
  TEname <- rlang::enquo(TE)
  tedata <- data %>% filter(fusionclass == "valid") %>% filter(TEtype == !!TEname) %>% 
    pivot_longer(cols = c(contains("Ctrl"), contains("MERVL")), names_to = "sample", values_to = "CPM") %>% 
    mutate(sample2=sample) %>% separate(sample2, sep = "_",into = c("KD","stage","rep")) %>% 
    group_by(gene5, strand,chr, start5, end5, X13,gene3, start3,end3,TEtype,fusionclass,KD,stage) %>% 
    summarise(CPM = sum(CPM)/2) %>% ungroup() %>% 
    group_by(gene5, strand,chr, start5, end5, X13,gene3, start3,end3,TEtype,fusionclass) %>% 
    mutate(total = sum(CPM)) %>% ungroup() %>% filter(total >=0.6)
  datalabel <- tedata %>% group_by(stage) %>% 
    summarise(result_stage = wilcox_test(.,CPM ~ KD, exact = TRUE)$stage,
              result_p = wilcox_test(.,CPM ~ KD, exact = TRUE)$p,
              count2 = n()/2) %>% ungroup() %>%  
    filter(stage == result_stage) %>% 
    mutate(result2 = paste("italic(p) ==", signif(result_p, digits = 2)),
           count2 = paste("italic(n) ==", count2),
           star = if_else(result_p > 1e-2, "NS",
                          if_else(result_p > 1e-4, "*",
                                  if_else(result_p > 2.2e-16, "**","***"))))
  tedata %>% 
    ggplot(aes(x = KD, y =log2(CPM+1)))+
    geom_violin(aes(fill = KD))+
    geom_boxplot(width =.1, fill = "black", outlier.colour = NA)+
    facet_wrap(~stage,nrow = 1)+
    stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2)+
    geom_text(x = 1.1, y = 5.3, data = datalabel, aes(label = result2), parse = TRUE, hjust = 0)+
    geom_text(x = 1.1, y = 5, data = datalabel, aes(label = count2), parse = TRUE, hjust = 0) +
    geom_text(x = 1.3, y = 5.6, data = datalabel, aes(label = star), hjust = 0,size = 3.5) +
    theme_minimal() +
    labs(y = "log2(CPM + 1)") +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank()) +
    ggsave(paste0(output,TE,"_box.pdf"), width =4, height =4)
}
cal_box_TE(chimeraall2, TE = "MERVL-int")
cal_box_TE(chimeraall2, TE = "MT2_Mm")
cal_box_TE(chimeraall2, TE = "MTA_Mm")
cal_box_TE(chimeraall2, TE = "MTC-int")
cal_box_TE(chimeraall2, TE = "RLTR45-int")
cal_box_TE(chimeraall2, TE = "IAPEz-int")



chimeraall2 %>% 
  filter(fusionclass == "valid") %>%
  write_csv(paste0(output,"all_chimera_list.csv"))











