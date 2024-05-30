
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/CLASH/vienna/")
suppressMessages(suppressWarnings(require(tidyverse)))

# count making!!!!!!!!!!!!!!!!!!!
ppi <- 300
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

outdir <- "/media/ariel/Isoseq/chimera_view/"

#######################################################################################################
# table making
#######################################################################################################
isoseq_count <- read_tsv("/media/ariel/Isoseq/process/polish_cluster_count.tsv")


count_computing <- function(samplename){
  isoform_matrix <- read_tsv(paste0(outdir, "polish_",samplename,"_selected.bed"),
                             col_names = c("gene","start","end","cluster_id","score","strand","dis1","dis2",
                                           "dis3","dis4","dis5","dis6")) %>% 
    select(-contains("dis"))
  
  isoform_matrix2 <- isoform_matrix %>% select(gene,cluster_id) %>% 
    left_join(isoseq_count %>% 
                pivot_wider(names_from = "sample",values_from = "count",values_fill = list(count = 0)) , by = "cluster_id") %>% 
    mutate(EGFP = siEGFP + siEGFP2, PIWI = siPIWI + siPIWI2) %>% select(-siEGFP,-siEGFP2, -siPIWI, -siPIWI2)
  isoform_matrix2 %>% 
    write_csv(paste0(outdir, "polish_",samplename,"_count.csv"))
}

count_computing(samplename = "Fuca_hap1")
count_computing(samplename = "CG7460_hap1")
count_computing(samplename = "Oat_hap2")
count_computing(samplename = "Shal_hap1")

# CAGE count calculation


CAGE_count <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "CAGE_Matrix.csv"))%>% arrange(sample) %>% 
  pivot_wider(names_from = "sample", values_from = "count", values_fill = list(count = 0)) %>% 
  mutate(EGFP = siEGFP_rep1 + siEGFP_rep2 + siEGFP_rep3, PIWI = siPIWI_rep1 + siPIWI_rep2 + siPIWI_rep3) %>% select(-contains("rep"))

CAGE_chimera <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/","CAGE_all_chimera_list.csv"))%>% arrange(sample) %>% 
  pivot_wider(names_from = "sample", values_from = "count", values_fill = list(count = 0)) %>% 
  mutate(EGFP = siEGFPrep1 + siEGFPrep2 + siEGFPrep3, PIWI = siPIWIrep1 + siPIWIrep2 + siPIWIrep3) %>% 
  mutate(fusiongroup = if_else(str_detect(chr_donorA, "chr"), "genome_TE", "TE_genome")) %>% 
  group_by(symbol, TEtype,fusiongroup) %>% summarise(EGFP = sum(EGFP), PIWI = sum(PIWI)) %>% ungroup()



TSS_chimera_count_cal <- function(genesymbol2, samplename, TEname) {
  CAGE_gene_count <- CAGE_count %>% 
    filter(genesymbol == genesymbol2) 
  if (dim(CAGE_gene_count)[1] != 0){
  CAGE_gene_count %>% write_csv(paste0(outdir, "TSS_count/TSS_gene_",samplename, "_",genesymbol2,"_count.csv"))
  }
  
  CAGE_chimera_count <- CAGE_chimera %>% 
    filter(symbol == genesymbol2) %>% filter(TEtype == TEname)
  if (dim(CAGE_chimera_count)[1] != 0){
  CAGE_chimera_count %>% write_csv(paste0(outdir, "TSS_count/TSS_chimera_",samplename, "_",genesymbol2,"_count.csv"))
  }
}


TSS_chimera_count_cal(genesymbol2 = "Fuca", samplename= "Fuca", TEname= "Dmel-17.6")
TSS_chimera_count_cal(genesymbol2 = "CG11714", samplename= "Fuca", TEname= "Dmel-17.6")
TSS_chimera_count_cal(genesymbol2 = "CG7460", samplename= "CG7460", TEname= "Dmel-mdg1")
TSS_chimera_count_cal(genesymbol2 = "Oat", samplename= "Oat", TEname= "Dmel-297")
TSS_chimera_count_cal(genesymbol2 = "CG9231", samplename= "Shal", TEname= "Dmel-gypsy")
TSS_chimera_count_cal(genesymbol2 = "Shal", samplename= "Shal", TEname= "Dmel-gypsy")

