
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gt)))


ppi <- 600
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

outdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/fig3/"

#######################################################################################################
# Fig3A, TE change large tables
#######################################################################################################
CAGE_table <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "CAGE_Matrix.csv")) %>% 
  select(genename, sample,count) %>% 
  group_by(genename, sample) %>% summarise(count = sum(count)) %>% ungroup() %>% arrange(genename,sample) %>% 
  pivot_wider(names_from = "sample",values_from = "count" , values_fill =list(count=0)) %>% 
  select(genename, siEGFP_rep1, siEGFP_rep2,siEGFP_rep3,siPIWI_rep1,siPIWI_rep2,siPIWI_rep3) %>% 
  mutate(total = rowSums(across(contains("rep"))),
         total2 = sum(total),
         total_CPM = 1000000*total/total2) %>% 
  filter(total_CPM >= 0.5) %>% 
  select(-total, -total2,-total_CPM) %>% 
  mutate(across(contains("rep"), ~(1000000 * .x/sum(.x)))) %>% 
  mutate(siEGFP = (siEGFP_rep1 + siEGFP_rep2 + siEGFP_rep3)/3,
         siPIWI = (siPIWI_rep1 + siPIWI_rep2 +siPIWI_rep3)/3,
         target = genename %>% str_replace("-", "|") %>% str_replace("_", "|"))  

Vienna_struct <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/", "CLASH_sup_file4.csv"),
                          col_types = cols(dG = col_double()))

Vienna_cluster <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig2/", "clash_selection_cluster.csv"), 
                           col_types = cols(dG = col_double()))

table4 <- read_tsv(paste0("/media/pericles/CLASH/ChIPseq/change/","change_edgeR_","H3K9me3","_","keio1","_","1000",".tsv"))

tx3_CLIP <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "CLIP_Matrix.csv")) %>% 
  filter(targettype =="ensemblTE") %>% separate(target, sep = "%", into = c("target", "dis1","dis2")) 
  
CLASH_struct <- Vienna_cluster %>%
  left_join(Vienna_struct %>% select(sample,readname,piRNA,struct,pilen),by = c("sample","readname","piRNA")) %>%
  filter(targettype =="ensemblTE") %>% separate(target, sep = "%", into = c("target", "dis1","dis2")) %>% 
  mutate(target = target %>% str_replace("element","-element") %>% str_replace("Beagle","-Beagle") %>% str_replace("TART", "TART-")) 
type_order <- c("perfect match","near-perfect match","imperfect match")

cov_Isoseq <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "Isoseq_Matrix.csv")) %>% 
  mutate(target = associated_gene %>% str_replace("-", "|") %>% str_replace("_", "|")) %>%
  mutate(across(c(fl_EGFP,fl_PIWI), ~(1000000 * .x/sum(.x)))) %>% 
  rename(siEGFP = fl_EGFP, siPIWI = fl_PIWI, genename = associated_gene) %>% select(genename,target,siEGFP,siPIWI)

wobble_all <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/","CLASH_wobble_annotation.csv"))


RNA_color_box <- function(data,RNAdata,RNAname,siEGFP,siPiwi,out,threshold = 0.5,width = 10, height = 10) {
  siEGFP <- rlang::enquo(siEGFP)
  siPiwi <- rlang::enquo(siPiwi)
  
  CLASH_struct2 <- data %>% 
    left_join(RNAdata,by = c("target")) %>% select(-contains("dis")) %>% 
    filter(genename %in% c("Dmel-mdg1","Dmel-412","Dmel-297","Dmel-springer","Dmel-gypsy","Dmel-Stalker","Dmel-Stalker2","Dmel-Stalker4",
                           "Dmel-copia","Dmel-3S18","Dmel-Tabor","Dmel-roo","Dmel-blood","Dmel-Quasiomodo","Dmel-flea","Dmel-Doc",
                           "Dmel-ZAM","Dmel-I-element","Dmel-17.6","Dmel-Juan","Dmel-Idefix","Dmel-F-element", "Dmel-invader4",
                           "Dmel-gypsy4","Dmel-mdg3","Dmel-gtwin","Dmel-diver","Dmel-gypsy12"))
  CLASH_struct3 <- CLASH_struct2 %>% mutate(picount= str_count(piRNAmatch,"\\(")) %>% 
    left_join(wobble_all %>% mutate(wobble = "wobble"), by = c("sample", "readname")) %>% mutate(wobble = if_else(is.na(wobble), "normal",wobble)) %>%  
    mutate(type = case_when(((!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)"))) & wobble =="normal") ~ "perfect match",
                            ((!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)"))) & wobble =="wobble") ~ "near-perfect match",
                            (str_detect(piRNAmatch,"^\\.\\(\\(\\(") & pilen ==(picount + 1)) &!(str_detect(targetmatch,"\\)\\.+\\)"))  ~ "near-perfect match",
                            TRUE ~ "imperfect match"))  %>% 
    group_by(genename, type) %>% count() %>% ungroup()
  g1 <-  CLASH_struct2 %>% # filter(!!siEGFP != 0 & !!siPiwi !=0) %>%
    filter(!!siEGFP >= threshold) %>%
    mutate(dif = log2(!!siPiwi/!!siEGFP)) %>%
    filter(targettype =="ensemblTE") %>% filter(is.finite(dG)) %>%
    ggplot(aes(x = dG,y = fct_reorder(genename,dif),color = dif))+
    geom_jitter(size = 0.5,alpha = 0.5, width = 0,height = 0.2)+
    geom_boxplot(color = "black", fill = "white", outlier.color = NA, alpha = 0)+
    theme_minimal()+ labs(x = "\u0394G (kcal/mol)", color = "log2FC")+
    theme(axis.title.y = element_blank(), legend.position = "bottom")+
    #scale_color_gradientn(colours = RColorBrewer::brewer.pal(11, "RdGy") %>% rev()+
    scale_color_gradient(low = "blue",high = "red",na.value = "gray74")
  TEfac <- CLASH_struct2 %>% # filter(!!siEGFP != 0 & !!siPiwi !=0) %>%
    filter(!!siEGFP >= threshold) %>% mutate(dif = log2(!!siPiwi/!!siEGFP)) %>%
    filter(targettype =="ensemblTE") %>% filter(is.finite(dG)) %>%
    select(genename,dif) %>% distinct()
  g2 <- TEfac %>% left_join(RNAdata %>% select(genename,!!siEGFP, !!siPiwi),by = "genename") %>%
    rename(`EGFP KD` = !!siEGFP,`Piwi KD` = !!siPiwi) %>%
    pivot_longer(-c(genename, dif), names_to = "type2", values_to = "CPM") %>%
    ggplot() +
    geom_tile(aes(x =type2, fill = CPM, y=fct_reorder(genename,dif)))+
    geom_text(aes(x =type2,y=fct_reorder(genename,dif),label = if_else(is.infinite(CPM),"NA",if_else(CPM==0,"NA",""))))+
    scale_fill_gradient(low = "white", high = "dodgerblue4",na.value = "gray74",trans = "log10",breaks = c(1,100,10000))+
    labs(x = RNAname,fill="")+ theme_minimal()+
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), legend.position = "bottom",
          strip.background=element_blank(), legend.key.width = unit(0.4,"cm"))
  g3 <- TEfac %>% mutate(type2 = "log2FC") %>%
    ggplot() +
    geom_tile(aes(x =type2, fill = dif, y=fct_reorder(genename,dif)))+
    geom_text(aes(x =type2, y=fct_reorder(genename,dif),label = if_else(is.infinite(dif),"NA",if_else(dif==0,"NA",""))))+
    labs(x = "",fill="")+ theme_minimal()+
    theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), legend.position = "bottom",
          strip.background=element_blank(), legend.key.width = unit(0.3,"cm"))+
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral") %>% rev(),
                         na.value = "gray74",breaks=c(-5,0,5),limits = c(-5,5),oob=scales::squish)
  g4 <- TEfac %>% mutate(target = genename %>% str_replace("-", "|") %>% str_replace("_", "|")) %>% 
    left_join(tx3_CLIP %>% select(target,PiwiCLIP_IN, PiwiCLIP_IP_rep1), by = "target") %>%
    rename(input = PiwiCLIP_IN, iCLIP = PiwiCLIP_IP_rep1) %>% select(-target) %>% 
    pivot_longer(-c(genename, dif), names_to = "type2", values_to = "CPM") %>%
    ggplot(aes(x =type2, fill = CPM, y=fct_reorder(genename,dif))) +
    geom_tile()+
    geom_text(aes(label = if_else(is.na(CPM),"NA",if_else(CPM==0,"NA",""))))+
    scale_fill_gradient(low = "white", high = "magenta4",trans = "log2", na.value = "gray74",breaks = c(2,256,65536))+
    labs(x = "iCLIP",fill="")+ theme_minimal()+
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), legend.position = "bottom",
          strip.background=element_blank(), legend.key.width = unit(0.4,"cm"))
  g5 <- TEfac %>%  mutate(target = genename %>% str_replace("-", "|") %>% str_replace("_", "|")) %>% 
    left_join(tx3_CLIP %>% select(target,PiwiCLIP_IN, PiwiCLIP_IP_rep1), by = "target") %>%
    mutate(type2 = "log2FC",dif2 = log2(PiwiCLIP_IP_rep1/PiwiCLIP_IN)) %>%
    ggplot(aes(x =type2, fill = dif2, y=fct_reorder(genename,dif))) +
    geom_tile()+
    geom_text(aes(label = if_else(is.na(dif2),"NA",if_else(is.infinite(dif2),"NA",""))))+
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral") %>% rev(),
                         na.value = "gray74",breaks=c(-5,0,5),limits = c(-5,5),oob=scales::squish)+
    labs(x = "",fill="")+ theme_minimal()+
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), legend.position = "bottom",
          strip.background=element_blank(), legend.key.width = unit(0.3,"cm"))
  g6 <- TEfac %>% 
    left_join(CLASH_struct3,by = "genename") %>%
    group_by(genename) %>% mutate(total = sum(n)) %>% ungroup() %>% 
    mutate(per = 100 * n/total) %>% 
    mutate(type = type %>% factor(levels = type_order)) %>% 
    ggplot(aes(y = fct_reorder(genename,dif), x = per, fill = type)) +
    geom_col(color ="black", position = "stack")+
    geom_text(x = 110,size = 2.5, aes(label = paste0("n=",total)))+ guides(fill = guide_legend(nrow = 3,byrow = TRUE))+
    xlab("piRNA pairing (%)") + theme_minimal() +scale_fill_manual(values = c("maroon1","royalblue1","gray65"))+
    scale_x_continuous(breaks = c(0,33.3,66.7,100), labels = c(0,33.3,66.7,100), limits = c(0,120))+
    theme(axis.title.y = element_blank(), legend.position = "bottom", axis.text.y = element_blank(),
          legend.title = element_blank(), panel.grid.minor.x = element_blank(),
          legend.key.size = unit(1, "cm"),legend.key.height = unit(0.05, 'cm'),legend.key.width = unit(0.5, 'cm'))
  graph <- cowplot::plot_grid(g1, g2, g3,g4,g5,g6, align = "h", nrow = 1,axis = "tb",
                              rel_widths = c(3/8,1/8,1/16,1/8,1/16,1/4))
  png(paste0(outdir,"RNA_change_",out,"_add.png"), width = width*ppi, height = height*ppi, res = ppi)
  print(graph)
  dev.off()
}
CLASH_struct %>% RNA_color_box(siEGFP=siEGFP,siPiwi=siPIWI,RNAdata = CAGE_table,RNAname = "CAGE",
                                    threshold = 0,out = "CAGE",width = 10, height = 5)



