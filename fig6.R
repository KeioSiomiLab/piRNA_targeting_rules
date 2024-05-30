
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gt)))


ppi <- 600
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

outdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/fig6/"

# wobble annotation

Vienna1_wobble <- read_tsv(paste0("/media/hermione/piRNAmodel/vienna/","All_model_vienna.tsv"))

wobble_pi <- Vienna1_wobble  %>% select(readname,  piRNAseq, piRNAmatch) %>% 
  mutate(piRNAmatch2 = piRNAmatch %>% str_replace_all("\\.","0") %>% str_replace_all("\\(","1") %>% str_split(pattern = ""),
         piRNAseq2 = piRNAseq %>% str_split(pattern = ""), index2 = row_number()) %>% 
  mutate(base = map(piRNAmatch2, function(x) {seq(1, length(x), by = 1)})) %>% 
  select(readname, index2, base, piRNAmatch2, piRNAseq2) %>% 
  unnest_legacy() %>% filter(piRNAmatch2 == "1") %>% 
  group_by(readname, index2) %>% mutate(base2 = row_number()) %>% ungroup() %>% select(-base, -piRNAmatch2)

wobble_tag <- Vienna1_wobble %>% select(readname, targetseq, targetmatch) %>% 
  mutate(targetmatch2 = targetmatch %>% stringi::stri_reverse() %>% str_replace_all("\\.","0") %>% str_replace_all("\\)","1") %>% str_split(pattern = ""),
         targetseq2 = targetseq %>% stringi::stri_reverse() %>% str_split(pattern = ""), index2 = row_number()) %>% 
  mutate(base = map(targetmatch2, function(x) {seq(1, length(x), by = 1)})) %>% 
  select(readname, index2, base, targetmatch2, targetseq2) %>% 
  unnest_legacy() %>% filter(targetmatch2 == "1") %>% 
  group_by(readname, index2) %>% mutate(base2 = row_number())%>% ungroup() %>% select(-base, -targetmatch2)

wobble_all <- full_join(wobble_pi, wobble_tag, by = c("readname", "index2","base2")) %>% 
  mutate(wobble = if_else((piRNAseq2 == "G" & targetseq2 == "T") | (piRNAseq2 == "T" & targetseq2 == "G"), "wobble","pair" )) %>% 
  filter(wobble == "wobble") %>% select(readname) %>% distinct()
wobble_all %>% write_csv(paste0(outdir,"model_wobble_annotation.csv"))



#######################################################################################################
# model, this is not used!
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
  mutate(WT_CPM = (siEGFP_rep1 + siEGFP_rep2 + siEGFP_rep3)/3,
         PIWI_CPM = (siPIWI_rep1 + siPIWI_rep2 + siPIWI_rep3)/3,
         FC = (siPIWI_rep1 + siPIWI_rep2 +siPIWI_rep3)/(siEGFP_rep1 + siEGFP_rep2 +siEGFP_rep3))
RNAtype <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_RNAtype.tsv")


CAGE_table2 <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "CAGE_edgeR.csv"))


remove_RNA <- c("ex","KCNQ","bru2","kel","hdc","CG46467","CG17716")
too_long_genes <- c("l(3)80Fg","Myo81F","Pzl","Parp","spok","Snap25","nvd","CG40470","AGO3","lovit","l(3)80Fl","l(3)80Fj","lt",
                    "Dbp80","DIP-lambda","CG17684","dpr21","CG41520","l(2)41Ab","rl","CG40006","Marf1","Maf1","Gprk1")
chr_other_genes <- read_tsv("/media/pericles/CLASH/database/anno/gene_sort.bed",
                            col_names = c("chr","start","end","genename","score","strand")) %>% 
  select(chr,genename) %>% filter(!(chr %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")))
# chr2L : "CG40006","Marf1","Gprk1","lt"
# chr2R : "DIP-lambda","CG17684","dpr21","CG41520","l(2)41Ab","rl","Maf1"
# chr3L : "l(3)80Fg","Snap25","nvd","CG40470","AGO3","lovit","l(3)80Fl","l(3)80Fj","Dbp80"
# chr3R : "Myo81F","Pzl","Parp","spok"
wobble_all <- read_csv(paste0(outdir,"model_wobble_annotation.csv"))
Vienna1_model <- read_tsv(paste0("/media/hermione/piRNAmodel/vienna/","All_model_vienna.tsv")) %>% 
  anti_join(chr_other_genes %>% rename(targetname = genename), by = "targetname") %>% 
  left_join(RNAtype %>% rename(RNAtype = type,targetname = name), by = "targetname") %>% 
  mutate(RNAtype = if_else(is.na(RNAtype), "TE",RNAtype),
         genesymbol = if_else(is.na(genesymbol), targetname,genesymbol),
       #  remover = if_else(genesymbol %in% remove_RNA, "remove","stay"),
         remover2 = if_else(genesymbol %in% too_long_genes, "remove","stay")) %>% 
  filter(RNAtype != "ncRNA") %>% # filter(remover != "remove") %>% 
  filter(remover2 != "remove") %>% filter(type3 != "rmst") %>% 
  mutate(picount= str_count(piRNAmatch,"\\(")) %>% 
  left_join(wobble_all %>% mutate(wobble = "wobble"), by = c("readname")) %>% mutate(wobble = if_else(is.na(wobble), "normal",wobble)) %>%  
  mutate(pistruct = case_when(((!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)"))) & wobble =="normal") ~ "perfect match",
                              ((!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)"))) & wobble =="wobble") ~ "functional",
                              (str_detect(piRNAmatch,"^\\.\\(\\(\\(") & pilen ==(picount + 1)) &!(str_detect(targetmatch,"\\)\\.+\\)"))  ~ "functional",
                              TRUE ~ "not functional"))  %>% 
  filter(pistruct != "not functional")

merge_bed <- read_tsv(paste0("/media/hermione/piRNAmodel/vienna/","All_model_vienna_merge.bed"),
                      col_names = c("targetname","bedstart","bedend"))
merge_bed2 <-  merge_bed %>% mutate(bindlen = bedend - bedstart) %>% 
  group_by(targetname) %>% summarise(bindlen = sum(bindlen)) %>% ungroup()

wobble_per <- Vienna1_model %>% filter(wobble == "wobble") %>% 
  group_by(targetname) %>% summarise(wobbleCPM = sum(CPM)) %>% ungroup()

Vienna2 <- Vienna1_model %>% filter(type3 !="rmst") %>%
  mutate(ki = exp((-dG * 4.184 * 1000)/(8.3145 * 299.15)),
         kiE  = ki * CPM) %>% 
  group_by(targetname, RNAtype,genesymbol) %>% 
  summarise(psi = sum(kiE),
            omega = sum (CPM),meanpilen = mean(pilen),
            meandG = mean(dG), pinum = n()) %>% ungroup() %>% 
  left_join(CAGE_table %>% rename(targetname = genename), by = "targetname") %>% drop_na() %>% 
  left_join(merge_bed2, by = "targetname") %>% left_join(wobble_per, by = "targetname") %>% 
  mutate(phi = WT_CPM * psi, chi = phi / omega, log2FC = log2(FC), log10CPM = log10(WT_CPM),
         wobbleCPM = if_else(is.na(wobbleCPM), 0, wobbleCPM)) %>% 
  mutate(change = if_else(log2FC>1, "change","stay"),
         tau = omega/WT_CPM,
         wobble_per = wobbleCPM/ omega,
         piRNA_density = omega / bindlen)
  

#######################################################################################################
# genic piRNA comparison
#######################################################################################################
# genic piRNA
fusion_target <- c("Oat","CG42565","CG15209","CHKov1","CG2233","NT5E-2","CG8303","TwdlL","Adgf-A","Irk1","CG7460","Cyp12e1")
fusion_target_selected <- c("Cyp6a9","kel","KCNQ","Rbp6","p130CAS","ogre","grp","HnRNP-K","Elk","pico","Ge-1")

RNAtype <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_RNAtype.tsv")
too_long_genes <- c("l(3)80Fg","Myo81F","Pzl","Parp","spok","Snap25","nvd","CG40470","AGO3","lovit","l(3)80Fl","l(3)80Fj","lt",
                    "Dbp80","DIP-lambda","CG17684","dpr21","CG41520","l(2)41Ab","rl","CG40006","Marf1","Maf1","Gprk1")
chr_other_genes <- read_tsv("/media/pericles/CLASH/database/anno/gene_sort.bed",
                      col_names = c("chr","start","end","genename","score","strand")) %>% 
  select(chr,genename) %>% filter(!(chr %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")))
  
genicpiRNA <- read_tsv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/piRNA/", "piRNA_CPM.tsv")) 
wobble_all <- read_csv(paste0(outdir,"model_wobble_annotation.csv"))
Vienna1 <- read_tsv(paste0("/media/hermione/piRNAmodel/vienna/","All_model_vienna.tsv")) %>% 
  anti_join(chr_other_genes %>% rename(targetname = genename), by = "targetname") %>% 
  left_join(RNAtype %>% rename(RNAtype = type,targetname = name), by = "targetname") %>% 
  mutate(RNAtype = if_else(is.na(RNAtype), "TE",RNAtype),
         genesymbol = if_else(is.na(genesymbol), targetname,genesymbol),
         remover2 = if_else(genesymbol %in% too_long_genes, "remove","stay")) %>% 
  filter(RNAtype != "ncRNA") %>%  filter(remover2 != "remove") %>% filter(type3 != "rmst") %>% 
  mutate(picount= str_count(piRNAmatch,"\\(")) %>% 
  left_join(wobble_all %>% mutate(wobble = "wobble"), by = c("readname")) %>% mutate(wobble = if_else(is.na(wobble), "normal",wobble)) %>%  
  mutate(pistruct = case_when(((!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)"))) & wobble =="normal") ~ "perfect match",
                              ((!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)"))) & wobble =="wobble") ~ "functional",
                              (str_detect(piRNAmatch,"^\\.\\(\\(\\(") & pilen ==(picount + 1)) &!(str_detect(targetmatch,"\\)\\.+\\)"))  ~ "functional",
                              TRUE ~ "not functional"))  %>% 
  filter(pistruct != "not functional")

genic_Vienna <- Vienna1 %>% 
  left_join(genicpiRNA %>% rename(piRNAname = name, pitype=type3,piRNAfrom=target) %>%
              select(piRNAname,pitype,piRNAfrom), by = "piRNAname") 
# stats of piRNA-target 
piRNA_class_order <- c("flam","20A","TE","genic","others")
RNA_class_order <- c("TE","mRNA","pseudogene")
genic_Vienna %>% mutate(pitype = pitype %>% factor(levels = piRNA_class_order),
                        RNAtype = RNAtype %>% factor(levels = RNA_class_order)) %>% 
  group_by(pitype,RNAtype) %>% count() %>% ungroup() %>% 
  ggplot(aes(x = pitype, y = RNAtype, fill = n)) + geom_tile()+geom_text(aes(label = n), size=2.5)+
  labs(x = "piRNA class",fill="", y = "target RNA type")+theme_minimal() +
  theme()+
  scale_fill_gradient(low = "white", high = "orange", na.value = "gray74",
                       breaks = c(0,10000,50000,100000), labels = c(0,10000,50000,100000)) +
  ggsave(paste0(outdir,"piRNA_RNA_target_type.png"), width =4.2, height =1.8, dpi = 600)+
  ggsave(paste0(outdir,"piRNA_RNA_target_type.pdf"), width =4.2, height =1.8)

genic_Vienna_mRNA_dist <- genic_Vienna %>% filter(RNAtype == "mRNA") %>% 
  group_by(targetname, genesymbol,pitype) %>% summarise(CPM = sum(CPM)) %>% ungroup() %>% 
  group_by(targetname) %>% mutate(CPM2 = sum(CPM)) %>% ungroup()
genic_Vienna_mRNA_list <- genic_Vienna_mRNA_dist %>% select(targetname, CPM2) %>%
  distinct() %>% arrange(desc(CPM2)) %>% head(50)

piRNA_class_order2 <- c("flam","TE","genic","others")
genic_Vienna_mRNA_list %>% select(targetname) %>% left_join(genic_Vienna_mRNA_dist, by = "targetname") %>% 
  mutate(pitype = pitype %>% factor(levels = piRNA_class_order2)) %>% 
  ggplot(aes(y = fct_reorder(genesymbol,CPM2), x = CPM, fill = pitype))+ geom_col()+
  theme_minimal()+ labs(x = "functional piRNA's total expression (CPM)")+
  theme(axis.title.y = element_blank(), legend.position = c(0.8,0.2),legend.title = element_blank())+
  ggsave(paste0(outdir,"piRNA_gene_expression.png"), width =4, height =6, dpi = 600)+
  ggsave(paste0(outdir,"piRNA_gene_expression.pdf"), width =4, height =6)

# stat of piRNA simulations

genic_Vienna %>% 
  group_by(pistruct) %>% count() %>%  write_csv(paste0(outdir, "stat_of_piRNA_simulations1.csv"))
genic_Vienna %>% 
  group_by(RNAtype) %>% count()%>%  write_csv(paste0(outdir, "stat_of_piRNA_simulations2.csv"))

# 

CAGE_table3 <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "CAGE_Matrix.csv")) %>% 
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
         siPIWI = (siPIWI_rep1 + siPIWI_rep2 +siPIWI_rep3)/3) %>%
  left_join(RNAtype %>% select(name, genesymbol) %>% rename(genename = name), by ="genename") %>% 
  mutate(genesymbol = if_else(is.na(genesymbol), genename, genesymbol))

piRNAexpress <- genic_Vienna %>% mutate(pitype = if_else(pitype %in% c("flam","TE","20A"), "TE",pitype)) %>%
  filter(pitype %in% c("TE","genic")) %>% select(genesymbol,pitype,CPM, pistruct) %>%
  left_join(CAGE_table3 %>% select(genesymbol,siEGFP,siPIWI), by = "genesymbol") %>% drop_na() %>% select(-siEGFP,-siPIWI) %>% 
  group_by(genesymbol,pitype) %>% summarise(CPM = sum(CPM),per = sum(pistruct == "perfect match")/n()) %>% ungroup() %>% 
  group_by(pitype) %>% slice_max(n=10, CPM, with_ties = FALSE) %>% ungroup()

CAGE_up <- enframe(fusion_target,value = "genesymbol") %>% mutate(type2 = "up")
CAGE_up_selected <- enframe(fusion_target_selected,value = "genesymbol") %>% mutate(type2 = "selected")

selectedTEgene <- bind_rows(CAGE_up,CAGE_up_selected) %>% select(-name) %>% 
  left_join(genic_Vienna, by = "genesymbol") %>% select(genesymbol,pitype,CPM,type2) %>% mutate(CPM = if_else(is.na(CPM), 0,CPM)) %>% 
  group_by(genesymbol,pitype,type2) %>% summarise(CPM = sum(CPM)) %>% ungroup()

table_ref <- piRNAexpress %>% mutate(type2 =if_else(pitype == "TE","TE related piRNAs\n(flam, 20A included)","genic piRNAs")) %>% 
  bind_rows(selectedTEgene %>% mutate(type2 =if_else(type2 == "up","Up regulated in CAGE\n(FDR top12)","Up regulated in CAGE\n(selected)"))) %>% 
  mutate(type2 = type2 %>% factor(levels = c("TE related piRNAs\n(flam, 20A included)","genic piRNAs",
                                             "Up regulated in CAGE\n(FDR top12)","Up regulated in CAGE\n(selected)"))) %>%
  rename(piCPM = CPM) %>% filter(type2 != "genic piRNAs")

table_ref %>% ggplot(aes(y = fct_reorder(genesymbol,piCPM), x = piCPM, fill = type2))+
  geom_col(color = "black")+ scale_x_sqrt(breaks = c(0,1,10,100,1000,10000,30000)) + guides(fill = "none")+
  facet_grid(type2 ~ .,  scales = "free_y",space = "free_y", switch = 'y')+
  theme_minimal()+ labs(y = "targeted gene",x = "functional piRNA's total expression (CPM)") +
  scale_fill_manual(values = c("aquamarine3","mediumpurple","maroon3","orange"))+ theme(strip.placement = "outside")+
  ggsave(paste0(outdir,"genic_expression.png"), width =4, height =8, dpi = 600)+
  ggsave(paste0(outdir,"genic_expression.pdf"), width =4, height =8)


piRNA_express_CAGE <- function(data = table_ref,RNAdata,siEGFP,siPiwi,out,width = 10, height = 10) {
  siEGFP <- rlang::enquo(siEGFP)
  siPiwi <- rlang::enquo(siPiwi)
g1 <- data  %>% ggplot(aes(y = fct_reorder(genesymbol,piCPM), x = piCPM, fill = type2))+
  geom_col(color = "black")+ scale_x_sqrt(breaks = c(0,100,1000,10000,30000),labels = c(0,100,1000,10000,30000)) +
  guides(fill = "none")+
  facet_wrap(~type2, ncol = 1, scales = "free_y",strip.position ="left")+
  theme_minimal()+ labs(y = "targeted gene",x = "piRNA expression (CPM)\n(perfect and near-perfect match)") +
  scale_fill_manual(values = c("aquamarine3","mediumpurple","maroon3","orange"))+
  theme(strip.placement = "outside", axis.text.x = element_text(angle = 40, hjust = 1), panel.grid.minor.x = element_blank())
g2 <- data  %>% ggplot(aes(y = fct_reorder(genesymbol,piCPM), x = 100*per) )+
  geom_col(color = "black",fill = "maroon1")+  
  facet_wrap(~type2, ncol = 1, scales = "free_y")+scale_x_continuous(breaks = c(0,50,100))+
  theme_minimal()+ labs(y = "",x = "perfect match (%)") +
  theme(strip.background=element_blank(),strip.text.x = element_blank(), axis.text.y = element_blank(),axis.title.y = element_blank())
g3 <- data  %>% select(genesymbol, type2,piCPM) %>% left_join(RNAdata %>% select(genesymbol,!!siEGFP, !!siPiwi),by = "genesymbol") %>%
  rename(`EGFP KD` = !!siEGFP,`Piwi KD` = !!siPiwi) %>%
  pivot_longer(-c(genesymbol, type2,piCPM), names_to = "type3", values_to = "CPM") %>%
  ggplot() +
  geom_tile(aes(x =type3, fill = CPM, y=fct_reorder(genesymbol,piCPM)))+
  geom_text(aes(x =type3,y=genesymbol,label = if_else(is.na(CPM),"NA",if_else(is.infinite(CPM),"NA",if_else(CPM==0, "NA","")))))+
  scale_fill_gradient(low = "white", high = "dodgerblue4",na.value = "gray74",trans = "log10",breaks = c(1,100,10000))+
  labs(x = "CAGE",fill="")+ facet_wrap(~type2, ncol = 1, scales = "free_y")+
  theme_minimal()+
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "bottom", strip.background=element_blank(),legend.key.width = unit(0.4,"cm"),strip.text.x = element_blank())
g4 <- data  %>% select(genesymbol, type2,piCPM) %>% left_join(RNAdata %>% select(genesymbol,!!siEGFP, !!siPiwi) ,by = "genesymbol") %>% 
  mutate(log2FC = log2(!!siPiwi/!!siEGFP),type3 = "log2FC") %>%
  ggplot() +
  geom_tile(aes(x =type3, fill = log2FC, y=fct_reorder(genesymbol,piCPM)))+
  geom_text(aes(x =type3, y=genesymbol,label = if_else(is.na(log2FC),"NA",if_else(is.infinite(log2FC),"NA",if_else(log2FC==0, "NA","")))))+
  labs(x = "",fill="")+theme_minimal()+facet_wrap(~type2, ncol = 1, scales = "free_y")+
  theme(axis.title.y = element_blank(),  axis.ticks.y = element_blank(), axis.text.y = element_blank(),strip.text.x = element_blank(),
        legend.position = "bottom", strip.background=element_blank(), legend.key.width = unit(0.3,"cm"))+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral") %>% rev(),
                       na.value = "gray74",breaks=c(-5,0,5),limits = c(-5,5),oob=scales::squish)
graph <- cowplot::plot_grid(g1, g2, g3,g4, align = "h", nrow = 1,axis = "tb",
                            rel_widths = c(1/2,1/6,2/9,1/9))
png(paste0(outdir,"Comparison_piRNA_simulation1.png"), width = width*ppi, height = height*ppi, res = ppi)
print(graph)
dev.off()
}
piRNA_express_CAGE(siEGFP=siEGFP,siPiwi=siPIWI,RNAdata = CAGE_table3,out = "CAGE",width = 8, height =6)


#######################################################################################################
# genic vs TEdrived
#######################################################################################################

All_fusion <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig4/","Insertion_chimera_table.csv")) %>% distinct()

table_ref2 <- piRNAexpress %>% mutate(type2 =if_else(pitype == "TE","TE related piRNAs\n(flam, 20A included)","genic piRNAs")) %>% filter(pitype!= "TE") %>%
  bind_rows(selectedTEgene %>% mutate(type2 =if_else(type2 == "up","Up regulated in CAGE\n(FDR top12)","Up regulated in CAGE\n(selected)"))) %>% 
  mutate(type2 = type2 %>% factor(levels = c("genic piRNAs",
                                             "Up regulated in CAGE\n(FDR top12)","Up regulated in CAGE\n(selected)"))) %>%
  rename(piCPM = CPM)%>% filter(type2 != "genic piRNAs")

All_fusion_map <- table_ref2 %>% select(genesymbol, type2,piCPM) %>% 
  left_join(All_fusion %>% select(genesymbol,TE, siEGFP, siPIWI,FC,strand1,direction), by = "genesymbol") %>% 
  group_by(genesymbol) %>% mutate(counter = n()) %>% ungroup() %>% mutate(FC = log2(siPIWI/siEGFP)) %>% 
  mutate(remove = case_when(is.na(TE)~ "stay",
                            genesymbol == "Glut4EF" & TE == "Dmel-412" ~ "stay",
                            genesymbol == "CG2233" & TE == "Dmel-gypsy" ~ "stay",
                            genesymbol == "kel" & TE == "Dmel-Quasimodo" ~ "stay",
                            genesymbol == "Oat" & TE == "Dmel-297" ~ "stay",
                            genesymbol == "Rbp6" & TE == "Dmel-297" ~ "stay",
                            genesymbol == "Rbp6" & TE == "Dmel-gypsy" ~ "stay",
                            counter == 1~ "stay",
                            TRUE ~ "remove"))  %>% 
  filter(remove == "stay") %>% distinct() %>% 
  mutate(genesymbol = if_else(genesymbol == "Rbp6" & TE == "Dmel-gypsy", "Rbp6#2",genesymbol),
         strandness = if_else(strand1 == direction, "sense","antisense"))


piRNA_fusion <- function(data = All_fusion_map,out,width = 10, height = 10) {
  g1 <- data  %>% ggplot(aes(y = fct_reorder(genesymbol,piCPM), x = remove))+
    geom_text(aes(label = TE),size = 3)+
    theme_minimal() + xlab("Neighboring TEs")+facet_wrap(~type2, ncol = 1, scales = "free_y")+
    theme(axis.title.y = element_blank(),panel.grid = element_blank(),axis.text.x = element_blank(),
          strip.text.x = element_blank(),strip.background=element_blank())
  g2 <- data  %>% ggplot(aes(y = fct_reorder(genesymbol,piCPM), x = remove))+
    geom_text(aes(label = strandness),size = 3)+
    theme_minimal() + xlab("strandedness to\nthe gene")+facet_wrap(~type2, ncol = 1, scales = "free_y")+
    theme(axis.title.y = element_blank(),panel.grid = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
          strip.text.x = element_blank(),strip.background=element_blank())
  g3 <- data  %>% rename(`EGFP KD` = siEGFP,`Piwi KD` = siPIWI) %>%pivot_longer(c(`EGFP KD`, `Piwi KD`), names_to = "type3", values_to = "CPM") %>%
    ggplot()+
    geom_tile(aes(x =type3, fill = CPM, y=fct_reorder(genesymbol,piCPM)))+
    geom_text(aes(x =type3,y=fct_reorder(genesymbol,piCPM),label = if_else(is.na(CPM),"NA",if_else(is.infinite(CPM),"NA",if_else(CPM==0, "NA","")))))+
    scale_fill_gradient(low = "white", high = "deeppink4",na.value = "gray74",trans = "log10",breaks = c(1,20,400))+
    labs(x = "chimera count",fill="")+ facet_wrap(~type2, ncol = 1, scales = "free_y")+
    theme_minimal()+
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          legend.position = "bottom", strip.background=element_blank(),legend.key.width = unit(0.4,"cm"),strip.text.x = element_blank())
  g4 <- data  %>% mutate(remove = "log2FC") %>% 
    ggplot() +
    geom_tile(aes(x =remove, fill = FC, y=fct_reorder(genesymbol,piCPM)))+
    geom_text(aes(x =remove, y=genesymbol,label = if_else(is.na(FC),"NA",if_else(is.infinite(FC),"NA",if_else(FC==0, "NA","")))))+
    labs(x = "",fill="")+theme_minimal()+facet_wrap(~type2, ncol = 1, scales = "free_y")+
    theme(axis.title.y = element_blank(),  axis.ticks.y = element_blank(), axis.text.y = element_blank(),strip.text.x = element_blank(),
          legend.position = "bottom", strip.background=element_blank(), legend.key.width = unit(0.3,"cm"))+
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral") %>% rev(),
                         na.value = "gray74",breaks=c(-5,0,5),limits = c(-5,5),oob=scales::squish)
  
  graph <- cowplot::plot_grid(g1, g2,g3,g4, align = "h", nrow = 1,axis = "tb",
                              rel_widths = c(4/9,1/6,2/9,1/6))
  png(paste0(outdir,"Comparison_piRNA_simulation2.png"), width = width*ppi, height = height*ppi, res = ppi)
  print(graph)
  dev.off()
}
piRNA_fusion(out = "CAGE",width = 6, height =5)

#######################################################################################################
# TE comparison
#######################################################################################################

piRNAexpressTE <- genic_Vienna %>% mutate(pitype = if_else(pitype %in% c("flam","TE","20A"), "TE",pitype)) %>%
  filter(pitype %in% c("TE") & RNAtype == "TE") %>% select(genesymbol,pitype,CPM, pistruct) %>%
  left_join(CAGE_table3 %>% select(genesymbol,siEGFP,siPIWI), by = "genesymbol") %>% drop_na() %>% select(-siEGFP,-siPIWI) %>% 
  group_by(genesymbol,pitype) %>% summarise(CPM = sum(CPM),per = sum(pistruct == "perfect match")/n()) %>% ungroup() 


table_ref3 <- piRNAexpressTE %>% rename(piCPM = CPM)


piRNA_express_CAGE_TE <- function(data = table_ref3,RNAdata,siEGFP,siPiwi,out,width = 10, height = 10) {
  siEGFP <- rlang::enquo(siEGFP)
  siPiwi <- rlang::enquo(siPiwi)
  g1 <- data  %>% ggplot(aes(y = fct_reorder(genesymbol,piCPM), x = piCPM))+
    geom_col(color = "black", fill = "aquamarine3")+ scale_x_sqrt(breaks = c(0,100,1000,10000,30000),labels = c(0,100,1000,10000,30000)) + guides(fill = "none")+
    theme_minimal()+ labs(y = "TE related piRNAs (flam, 20A included)",x = "piRNA expression (CPM)\n(perfect and near-perfect match)") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid.minor.x = element_blank(), axis.title.y = element_blank())
  g2 <- data  %>% select(genesymbol, piCPM) %>% left_join(RNAdata %>% select(genesymbol,!!siEGFP, !!siPiwi),by = "genesymbol") %>%
    rename(`EGFP KD` = !!siEGFP,`Piwi KD` = !!siPiwi) %>% mutate(pirate = piCPM / `EGFP KD`) %>% 
    ggplot(aes(y = fct_reorder(genesymbol,piCPM), x = pirate))+
    geom_col(color = "black", fill = "darkslateblue")+ scale_x_log10() + guides(fill = "none")+
    theme_minimal()+ labs(y = "",x = "piRNA/CAGE (CPM)") +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.grid.minor.x = element_blank())
  g3 <- data  %>% ggplot(aes(y = fct_reorder(genesymbol,piCPM), x = 100*per) )+
    geom_col(color = "black",fill = "maroon1")+scale_x_continuous(breaks = c(0,50,100))+
    theme_minimal()+ labs(y = "",x = "perfect match (%)") +
    theme(axis.text.y = element_blank(),axis.title.y = element_blank())
  g4 <- data  %>% select(genesymbol, piCPM) %>% left_join(RNAdata %>% select(genesymbol,!!siEGFP, !!siPiwi),by = "genesymbol") %>%
    rename(`EGFP KD` = !!siEGFP,`Piwi KD` = !!siPiwi) %>%
    pivot_longer(-c(genesymbol, piCPM), names_to = "type3", values_to = "CPM") %>%
    ggplot() +
    geom_tile(aes(x =type3, fill = CPM, y=fct_reorder(genesymbol,piCPM)))+
    geom_text(aes(x =type3,y=genesymbol,label = if_else(is.na(CPM),"NA",if_else(is.infinite(CPM),"NA",if_else(CPM==0, "NA","")))))+
    scale_fill_gradient(low = "white", high = "dodgerblue4",na.value = "gray74",trans = "log10",breaks = c(1,100,10000))+
    labs(x = "CAGE",fill="")+ theme_minimal()+
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          legend.position = "bottom", legend.key.width = unit(0.4,"cm"))
  g5 <- data  %>% select(genesymbol,piCPM) %>% left_join(RNAdata %>% select(genesymbol,!!siEGFP, !!siPiwi) ,by = "genesymbol") %>% 
    mutate(log2FC = log2(!!siPiwi/!!siEGFP),type3 = "log2FC") %>%
    ggplot() +
    geom_tile(aes(x =type3, fill = log2FC, y=fct_reorder(genesymbol,piCPM)))+
    geom_text(aes(x =type3, y=genesymbol,label = if_else(is.na(log2FC),"NA",if_else(is.infinite(log2FC),"NA",if_else(log2FC==0, "NA","")))))+
    labs(x = "",fill="")+theme_minimal()+
    theme(axis.title.y = element_blank(),  axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          legend.position = "bottom", legend.key.width = unit(0.3,"cm"))+
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral") %>% rev(),
                         na.value = "gray74",breaks=c(-5,0,5),limits = c(-5,5),oob=scales::squish)
  graph <- cowplot::plot_grid(g1, g2, g3,g4,g5, align = "h", nrow = 1,axis = "tb",
                              rel_widths = c(3/8,1/4,1/8,1/6,1/12))
  png(paste0(outdir,"Comparison_piRNA_simulation3.png"), width = width*ppi, height = height*ppi, res = ppi)
  print(graph)
  dev.off()
}
piRNA_express_CAGE_TE(siEGFP=siEGFP,siPiwi=siPIWI,RNAdata = CAGE_table3,out = "CAGE",width = 10, height =8)


#######################################################################################################
# CLASH result comparison
#######################################################################################################

Vienna_select <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/","CLASH_sup_file2.csv"),
                          col_types = cols(dG = col_double())) %>% filter(piRNAtype != "rmstpiRNA") %>% filter(targettype %in% c("mRNA","ensemblTE")) %>% 
  select(sample, readname, piRNA,target,targettype, symbol,dG) %>% 
  separate(target, sep = "%", into = c("target", "dis1","dis2")) %>% select(-contains("dis")) %>% 
  mutate(target = target %>% str_replace("element","-element") %>% str_replace("Beagle","-Beagle") %>% str_replace("TART", "TART-") %>% 
           str_replace("\\|","-")) %>% distinct()

CLASH_result <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/","CLASH_sup_file1.csv"), 
                          col_types = cols(dG = col_double())) %>% 
  filter(targettype%in% c("mRNA","ensemblTE")) %>% 
  mutate(picount= str_count(piRNAmatch,"\\(")) %>% 
  mutate(type = case_when(!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)")) ~ "perfect match",
   (str_detect(piRNAmatch,"^\\.\\(\\(\\(") & pilen ==(picount + 1)) &!(str_detect(targetmatch,"\\)\\.+\\)"))  ~ "functional",
                              TRUE ~ "not functional"))  %>% 
  select(sample,readname,piRNA,targettype,piRNAtype,type) %>%  
  left_join(Vienna_select, by = c("sample","readname","piRNA","targettype")) %>% drop_na() %>% distinct()

CLASH_result2 <- CLASH_result %>% 
  filter(type != "not functional") %>% 
  anti_join(chr_other_genes %>% rename(target = genename), by = "target") %>% 
  mutate(remover2 = if_else(symbol %in% too_long_genes, "remove","stay")) %>% filter(remover2 != "remove") %>% 
  left_join(genicpiRNA %>% rename(piRNA = name, CPM = CPM_avg) %>% select(piRNA,CPM), by = "piRNA") %>% 
  mutate(symbol = if_else(str_detect(target,"FBgn"),symbol,target))

CLASH_result3 <- CLASH_result2 %>%
  mutate(CLASHcount = readname %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer()) %>% 
  group_by(symbol) %>% summarise(CLASHcount  = sum(CLASHcount)) %>% ungroup()

Vienna_all <- Vienna1 %>% 
  group_by(genesymbol) %>% summarise(piCPM_sim = sum(CPM)) %>% ungroup() %>% 
  left_join(CAGE_table3 %>% select(genesymbol,siEGFP, siPIWI),by = "genesymbol") %>%
  rename(`siEGFP` = siEGFP,`siPiwi` = siPIWI) %>% 
  mutate(pirate = piCPM_sim / siEGFP,
         Teyes = if_else(genesymbol %in% paste0("Dmel-",dependent_TE), "Piwi dependent TE","Piwi independent"))


#######################################################################################################
# ROC results
#######################################################################################################

dependent_TE <- c("297","412","gypsy","gtwin","blood","17.6","mdg1","Tabor","I-element",
                  "HMS-Beagle2","rover","ZAM","Doc","flea","Quasimodo", "Stalker",
                  "Stalker2","Stalker4","Idefix","mdg3","F-element","gypsy5","Juan")

ROC_cal <- function(data, sorter,samplename) {
  sorter2 <- rlang::enquo(sorter)
  model_pi <- data %>% 
    arrange(desc(!!sorter2)) %>% 
    mutate(order = row_number() %>% as.integer(),
           tmp = if_else(Teyes == "Piwi depenedent TE",1,0),
           total_posi = sum(tmp),
           tmp2 = if_else(Teyes == "Piwi depenedent TE",0,1),
           total_nega = sum(tmp2),
           counter_posi = cumsum(tmp)/total_posi,
           counter_nega = cumsum(tmp2)/total_nega,
           sample = samplename) 
  return(model_pi)
}

Vienna3 <- Vienna2 %>% 
  mutate(Teyes = if_else(genesymbol %in% paste0("Dmel-",dependent_TE), "Piwi dependent TE","Piwi independent TE")) %>% 
  select(targetname, RNAtype, genesymbol,Teyes,psi,omega,phi,chi,tau,WT_CPM,PIWI_CPM, 
         meandG,pinum,meanpilen,bindlen, wobble_per, piRNA_density)
model_omega <- ROC_cal(data = Vienna3 ,sorter = omega,samplename = "omega")
model_tau <- ROC_cal(data = Vienna3 ,sorter = tau,samplename = "tau")

bind_rows(model_tau %>% mutate(samplename = "piRNA/CAGE"),
          model_omega %>% mutate(samplename = "piRNA expression")) %>%
  ggplot(aes(x = counter_nega, y = counter_posi, color = samplename))+
  geom_line(alpha = 0.8, size=1)+ geom_abline(aes(intercept = 0, slope = 1), color = "gray40",linetype = "dashed") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1))+scale_x_continuous(limits = c(0, 1))+
  labs(x = "False negative ratio",y = "True positive ratio") +
  theme(legend.title = element_blank(),legend.position = c(0.6,0.2),panel.grid.minor = element_blank()) +
  ggsave(paste0(outdir,"model_ROC_3.png"), width =2.5, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"model_ROC_3.pdf"), width =2.5, height =2.5)

AUC_cal <- function(data) {
  model_area <- data %>% filter(tmp2 == 1) %>% 
    mutate(area_diff = counter_posi / total_nega) %>% group_by(sample) %>%
    summarise(AUC = 100 * sum(area_diff)) %>% ungroup()
  return(model_area)
}
AUC_omega <- AUC_cal(data = model_omega)
AUC_tau <- AUC_cal(data = model_tau)

AUC_table <- bind_rows(AUC_omega%>% mutate(samplename = "piRNA expression"),
          AUC_tau %>% mutate(samplename = "piRNA/mRNA rate")) %>% 
  write_csv(paste0(outdir,"AUC_table.csv"))
  



Vienna3 %>% 
  ggplot(aes(x = omega, y =tau, color = Teyes))+
  geom_point(alpha = 0.8)+scale_x_log10()+scale_y_log10() + theme_minimal()+
  ggrepel::geom_text_repel(data = subset(Vienna3 %>% filter(genesymbol %in% c("Dmel-flea", "Dmel-Juan","Dmel-blood","Dmel-mdg3","Dmel-Doc"))),
                           aes(label =genesymbol), size = 2.7, color = "magenta4",  point.padding = 0.1,max.overlaps =Inf,nudge_x = -0.5, nudge_y =0.5,
                           box.padding = unit(1.2, "lines"), min.segment.length = 0)+
  ggrepel::geom_text_repel(data = subset(Vienna3 %>% filter(genesymbol %in% c("Dmel-springer","Dmel-copia","Dmel-gypsy4","Dmel-roo"))),
                           aes(label =genesymbol), size = 2.7, color = "black",  point.padding = 0.1,max.overlaps =Inf,nudge_x = 0.5, nudge_y =-0.5,
                           box.padding = unit(1.2, "lines"), min.segment.length = 0)+
  labs(x = "piRNA expression (CPM)\n(perfect and near-perfect match)",y = "piRNA/CAGE") +
  theme(legend.title = element_blank(),legend.position = c(0.25,0.9))+
  ggsave(paste0(outdir,"TE_piCPM_rate2.png"), width =3.5, height =3.5, dpi = 300)+
  ggsave(paste0(outdir,"TE_piCPM_rate2.pdf"), width =3.5, height =3.5)

#######################################################################################################
# cluster piRNA
#######################################################################################################

Vienna_cluster3 <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig2/","RNA_cluster_class.csv"),
            col_types = cols(dG = col_double()))

cluster_piRNA <- Vienna_cluster3 %>% 
  mutate(target = target %>% str_replace("\\|","-")) %>% 
  left_join(Vienna1 %>% group_by(targetname, RNAtype, genesymbol) %>% 
              summarise(CPM = sum(CPM)) %>% ungroup() %>% rename(target = targetname),by = "target")

cluster_piRNA %>% filter(!(is.na(CPM))) %>% filter(RNAtype != "pseudogene") %>% filter(cluster!="non") %>% 
  ggplot(aes(x = CPM, y = RNAtype,color=RNAtype)) +
  geom_jitter(size = 1.2,alpha = 0.8, width = 0,height = 0.2)+
  geom_boxplot(color = "black", fill = "white", outlier.color = NA, alpha = 0)+
  scale_x_log10()+
  facet_wrap(~ cluster, ncol = 1, strip.position="right")+
  theme_minimal()+ xlab("functional piRNA's expression level (CPM)")+ 
  theme(axis.title.y = element_blank())+ guides(color = "none")+
  scale_color_manual(values = c("goldenrod3","indianred1")) +
  ggsave(paste0(outdir,"cluster_piRNA_express.png"), width =3.5, height =2.8, dpi = 600)+
  ggsave(paste0(outdir,"cluster_piRNA_express.pdf"), width =3.5, height =2.8 )


#######################################################################################################
# GRO-seq cluster analysis
#######################################################################################################
#GRO-seq analysis
RNAtype <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_RNAtype.tsv")
too_long_genes <- c("l(3)80Fg","Myo81F","Pzl","Parp","spok","Snap25","nvd","CG40470","AGO3","lovit","l(3)80Fl","l(3)80Fj","lt",
                    "Dbp80","DIP-lambda","CG17684","dpr21","CG41520","l(2)41Ab","rl","CG40006","Marf1","Maf1","Gprk1")
GRO_gene <- read_tsv("/media/pericles/CLASH/RNAseq/GROcount/SRR609665_featurecount.txt",
                     col_names = c("name","len","exp"),skip = 2)
GRO_TE <- read_tsv("/media/pericles/CLASH/RNAseq/GROcount/SRR609665_TE.txt",
                   col_names = c("name","len","exp","dis")) %>% filter(len != 0) %>% select(-dis)
GRO_matrix <- bind_rows(GRO_gene,GRO_TE) %>% 
  mutate(CPM = 1000000*exp/sum(exp)) %>% 
  left_join(RNAtype,by = "name") %>%
  mutate(genesymbol = if_else(is.na(genesymbol), name,genesymbol),
         type = if_else(is.na(type), "TE",type)) %>% 
  filter(!(genesymbol %in% too_long_genes))

# CLASH cluster

Vienna_cluster <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig2/", "clash_selection_cluster.csv"))
Vienna_sle <-read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/", "CLASH_sup_file1.csv"),
                      col_types = cols(dG = col_double()))
test1 <- Vienna_sle %>% left_join(Vienna_cluster %>% select(readname, sample,targetname, targettype, target,cluster), 
                                  by = c("readname", "sample","targettype")) %>% 
  filter(cluster =="cluster3") %>% filter(targettype %in% c("mRNA")) %>% 
  group_by(sample, readname) %>% mutate(counter = 1/n()) %>% ungroup() %>% 
  group_by(target) %>% summarise(count = sum(counter)) %>% ungroup()

# figure
GRO_matrix %>% filter(type =="mRNA") %>% filter(CPM!=0) %>% 
  left_join(test1 %>% rename(name = target), by = "name") %>% drop_na() %>% 
  ggplot(aes(x = count, y= CPM)) +
  geom_point(alpha = 0.4)+geom_density_2d(alpha = 0.7, color = "pink")+
  geom_smooth(method = "lm")+
  scale_y_log10(breaks = c(0.1,1,10,100,1000),labels = c(0.1,1,10,100,1000))  +scale_x_log10()+
  theme_minimal()+ labs(y = "GRO-seq(CPM)",x = "CLASH count")+
  ggsave(paste0(outdir,"GRO_scatter_cluster3.png"), width =3.5, height =3.5, dpi = 300) +
  ggsave(paste0(outdir,"GRO_scatter_cluster3.pdf"), width =3.5, height =3.5)

GRO_matrix %>% filter(type =="mRNA") %>% filter(CPM!=0) %>% 
  left_join(test1 %>% rename(name = target), by = "name") %>% mutate(type = if_else(is.na(count), "unbound","bound")) %>% 
  ggplot(aes(x = CPM, fill= type))+
  geom_density(alpha = 0.5)+scale_x_log10(breaks = c(0.1,1,10,100,1000),labels = c(0.1,1,10,100,1000))+
  theme_minimal()+
  theme(legend.position = c(0.2,0.9),legend.title = element_blank())+labs(x = "GRO-seq(CPM)")+
  ggsave(paste0(outdir,"GRO_scatter_cluster3_distri.png"), width =3, height =3, dpi = 300) +
  ggsave(paste0(outdir,"GRO_scatter_cluster3_distri.pdf"), width =3, height =3)

referenceRNAorder <- c("TE","mRNA","ncRNA","rRNA","snoRNA","snRNA","tRNA","pseudogene")
GRO_matrix %>%  mutate(totalCPM = sum(CPM)) %>% filter(!(type%in% c("miRNA","pre_miRNA"))) %>% 
  group_by(type) %>% summarise(totalper = 100*sum(CPM)/totalCPM) %>% ungroup() %>% distinct() %>% 
  mutate(type = type %>% factor(levels = referenceRNAorder)) %>% 
  ggplot(aes(x="", y=totalper, fill=type)) +
  geom_bar(stat="identity", width=1, color="white", linewidth=0.1) +
  coord_polar("y", start=0, direction = 1) +
  theme_void() +
  guides(fill= guide_legend(title = NULL))+
  geom_text(aes(label = if_else(totalper > 3,paste0(round(totalper,1), "%"),"")), position = position_stack(vjust = 0.5))+
  ggsave(paste0(outdir,"GRO_pie.png"), width =4, height =3, dpi = 300) +
  ggsave(paste0(outdir,"GRO_pie.pdf"), width =4, height =3)

