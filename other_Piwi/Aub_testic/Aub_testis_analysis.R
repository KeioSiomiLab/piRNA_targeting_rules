

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gt)))
library(GGally)
library(ggthemes)
library(edgeR)
library(ggrepel)

ppi <- 600
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

cor_func <- function(data,mapping,method,symbol){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x,y, method = method,use = 'complete.obs')
  colFn <- colorRampPalette(c("brown1","white","dodgerblue"),interpolate = "spline")
  fill <- colFn(100)[findInterval(corr,seq(-1,1,length = 100))]
  ggally_text(label = paste(symbol,as.character(round(corr, 2))),
              mapping = aes(),
              xP = 0.5,yP = 0.5,color = "black", size = 2.5)+
    theme(panel.background = element_rect(fill = fill))
}
theme_set(theme_minimal())

RNAtype <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_RNAtype.tsv")


outdir <- "/media/hermione/otherdata_compare/Aub_testis/figures/"
readdir <- "/media/hermione/otherdata_compare/Aub_testis/"
#######################################################################################################
# RNA-seq analysis
#######################################################################################################
#mRNA
mRNA_dir <- dir(paste0(readdir, "featurecount/"), pattern = "\\_featurecount.txt$", full.names = TRUE)
mRNA_files <- mRNA_dir %>%
  str_remove(paste0(readdir, "featurecount/", "/")) %>% str_remove("\\_featurecount.txt")

df <- vector("list", length(mRNA_dir))
for (i in seq_along(mRNA_dir)) {
  df[[i]] <- read_tsv(mRNA_dir[[i]], skip = 2, col_names = c("gene_name", "length","count"))
  df[[i]][["type"]] <- mRNA_files[[i]]
}

mRNA_count <- bind_rows(df) %>% 
  select(gene_name,type,count) 
#TE
TE_dir <- dir(paste0(readdir, "TE/"), pattern = "\\_TEcount.txt$", full.names = TRUE)
TE_files <- TE_dir %>%
  str_remove(paste0(readdir, "TE/", "/")) %>% str_remove("\\_TEcount.txt")

df2 <- vector("list", length(TE_dir))
for (i in seq_along(TE_dir)) {
  df2[[i]] <- read_tsv(TE_dir[[i]], col_names = c("gene_name", "length","count","dis")) %>% 
    filter(gene_name != "*")
  df2[[i]][["type"]] <- TE_files[[i]]
}

TE_count <- bind_rows(df2) %>% 
  select(gene_name,type,count) 

read_info_RNA <- read_tsv(paste0(readdir,"read_info.tsv")) %>% 
  mutate(type = file %>% str_remove("_trim.fastq") %>% str_remove("_trimed3.fastq") %>%
           str_remove("Discas/") %>% str_remove("_1") %>% str_remove("_2") %>% str_remove("\\.gz"),
         num_seqs = num_seqs %>% str_remove_all(",") %>% as.double()) %>%
  select(type, num_seqs) %>% distinct()

All_count <- bind_rows(TE_count %>% mutate(RNAtype = "TE"), mRNA_count %>% mutate(RNAtype = "mRNA")) %>% 
  left_join(read_info_RNA, by = "type") %>% 
  mutate(CPM = count * 1000000/num_seqs) %>% select(-count, -num_seqs) %>%
  pivot_wider(names_from = type, values_from = CPM, values_fill =list(count=0))
# ggally
png(paste0(outdir,"cor_each_sample.png"), width = 5*ppi, height = 5*ppi, res = ppi)
pm <- ggpairs(data = All_count %>% select(-contains("gene_name")) %>%
                mutate_at(vars(contains("SRR")),list(~log10(.+1))),aes(color = RNAtype),
              xlab = "log10(target count + 1) ",ylab = "log10(target count + 1)",
              lower = list(continuous = wrap("points", alpha = 0.7, size = 0.2)),
              upper = list(continuous = wrap(cor_func, method = "pearson", symbol = expression('p ='))))+
  theme(text = element_text(size = 6), panel.grid.minor = element_blank())
pm2 <- pm
for(i in 2:pm$nrow) {
  for(j in 1:(i-1)) {
    pm2[i,j] <- pm[i,j] +
      scale_x_continuous(limits = c(0, 3.45)) +
      scale_y_continuous(limits = c(0, 3.45))}}
pm2
dev.off()

#edgeR

All_count %>% write_csv(paste0(outdir, "RNA_Matrix.csv"))
#edit here#####################################################################################################
All_edgeR <- bind_rows(TE_count %>% mutate(RNAtype = "TE"), mRNA_count %>% mutate(RNAtype = "mRNA")) %>% 
  select(gene_name, type,count) %>% pivot_wider(names_from = type, values_from = count, values_fill =list(count=0)) %>% 
  mutate(ctrl = SRR12216403, KD = (SRR12216401 + SRR12216402)/2) %>% select(-contains("SRR")) %>% 
  column_to_rownames(var = "gene_name") %>% 
  mutate(total = rowSums(across(everything())),
         total2 = sum(total),
         total_CPM = 1000000*total/total2) %>% 
  filter(total_CPM >= 0.5) %>% 
  select(-total, -total2,-total_CPM)

count <- All_edgeR
# multiple replicate
#group <- factor(c("ctrl","ctrl","KD","KD"))
#d <- DGEList(counts = count, group = group)
#d2 <- calcNormFactors(d, method = "RLE")
#design <- model.matrix( ~ 0 + group, data = d)
#d3 <- estimateGLMCommonDisp(d2, design)
#result <- exactTest(d3)
# single replicate
group <- factor(c("ctrl","KD"))
d <- DGEList(counts = count, group = group)
d2 <- calcNormFactors(d, method = "RLE")
d3 <- estimateGLMCommonDisp(d2, method = "deviance", robust = T, subset = NULL)
result <- exactTest(d3)

###############################################################################################################

table <- as.data.frame(topTags(result, n = nrow(count))) %>%
  rownames_to_column(var = "name") %>% as_tibble() %>%
  mutate(fdr2 = ifelse(FDR < 0.0001, "change", "stay"))

table2 <- table %>% left_join(All_count %>% rename(name = gene_name) %>% select(name, RNAtype), by = "name") %>% 
  left_join(count %>% rownames_to_column(var = "name"), by = "name") %>% 
  left_join(RNAtype, by = "name") %>% 
  mutate(type = if_else(is.na(type), "TE",type),
         genesymbol = if_else(is.na(genesymbol), name,genesymbol))
#edit here#####################################################################################################
changed_gene_list <- c("CG12717", "Ste:CG33237","Ste:CG33246","Ste:CG33247","Ste:CG33239","Ste:CG33242","Ste:CG33244","Ste:CG33245","Ste:CG33240","Ste:CG33241",
                       "Dmel-hopper2","Dmel-mdg3","Dsim-ninja","Dmel-aurora-element","Dmel-invader3","Dmel-invader2","Dmel-diver",
                       "Dmel-F-element","Dmel-X-element","Dmel-roo","Dmel-Doc2-element","Dmel-Xanthias","Dmel-Quasimodo",
                       "Dmel-invader6","Dmel-Circe","Dmel-412","Dmel-transib2","Dmel-Ivk","Dmel-Doc","Dmel-G6")
###############################################################################################################
change_order <- c("negative","positive","target RNA")
table3 <- table2 %>%
  mutate(ctrl_CPM = 1000000 * ctrl/ sum(ctrl),
         KD_CPM = 1000000 * KD/ sum(KD),
         per =  log2((KD_CPM+1)/(ctrl_CPM+1)),
         type = if_else(fdr2=="change", "positive", "negative"))  %>% 
  mutate(type = case_when(genesymbol %in% changed_gene_list ~ "target RNA",
                          TRUE ~ type)%>% factor(levels = change_order)) 

table3 %>%  
  ggplot(aes(x = log10(ctrl_CPM + 1), y = log10(KD_CPM + 1), color = type))+
  geom_point(size = 1.4,alpha = 0.6)+ scale_color_manual(values = c( "#808080","#dc143c","blue")) +
#  geom_text_repel(data = subset(table3, type == "target RNA"),
#                  aes(label = genesymbol ), color = "blue4", size = 2,max.overlaps = Inf,force=8, max.time = 1,
#                  segment.color = "gray90", fontface = "bold", segment.size = 0.5)+
  labs(x = "log10(control (CPM) + 1)",y = "log10(mutant (CPM) + 1)")+
  theme_minimal()+ guides(color = "none")+
  ggsave(paste0(outdir,"RNA_scatter.png"), width =3.5, height =3.5, dpi = 600)
table3 %>% write_csv(paste0(outdir, "RNA_edgeR.csv"))
#######################################################################################################
# piRNA analysis
#######################################################################################################
# wobble annotation

Vienna1_wobble <- read_tsv(paste0(readdir,"vienna/","All_model_vienna.tsv"))

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
# model
#######################################################################################################
remove_RNA <- c("ex","KCNQ","bru2","kel","hdc","CG46467","CG17716")
too_long_genes <- c("l(3)80Fg","Myo81F","Pzl","Parp","spok","Snap25","nvd","CG40470","AGO3","lovit","l(3)80Fl","l(3)80Fj","lt",
                    "Dbp80","DIP-lambda","CG17684","dpr21","CG41520","l(2)41Ab","rl","CG40006","Marf1","Maf1","Gprk1")
chr_other_genes <- read_tsv("/media/pericles/CLASH/database/anno/gene_sort.bed",
                            col_names = c("chr","start","end","genename","score","strand")) %>% 
  select(chr,genename) %>% filter(!(chr %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")))

wobble_all <- read_csv(paste0(outdir,"model_wobble_annotation.csv"))


gene_exon <- read_tsv(paste0(readdir ,"rnaplex/","gene_remove1_overlap.tsv"),
                        col_names = c("dis1","dis2","dis3","readname","targetname","dis4","dis5","TE","annotype")) %>% 
  filter(!(str_detect(annotype,"TE"))) %>% filter(annotype != "intron") %>% select(-contains("dis"))


Vienna1_model <- read_tsv(paste0(readdir,"vienna/","All_model_vienna.tsv")) %>% 
  left_join(gene_exon %>% select(readname, targetname,annotype) %>% distinct(), by = c("readname", "targetname")) %>% 
  filter(!(pairtype == "gene" & is.na(annotype))) %>% select(-annotype) %>% 
  anti_join(chr_other_genes %>% rename(targetname = genename), by = "targetname") %>% 
  left_join(RNAtype %>% rename(RNAtype = type,targetname = name), by = "targetname") %>% 
  mutate(RNAtype = if_else(is.na(RNAtype), "TE",RNAtype),
         genesymbol = if_else(is.na(genesymbol), targetname,genesymbol),
         #  remover = if_else(genesymbol %in% remove_RNA, "remove","stay"),
         remover2 = if_else(genesymbol %in% too_long_genes, "remove","stay")) %>% 
  filter(RNAtype != "ncRNA") %>% # filter(remover != "remove") %>% 
  filter(remover2 != "remove") %>% 
  mutate(picount= str_count(piRNAmatch,"\\(")) %>% 
  left_join(wobble_all %>% mutate(wobble = "wobble"), by = c("readname")) %>% mutate(wobble = if_else(is.na(wobble), "normal",wobble)) %>% 
  mutate(target_rev = targetmatch %>% str_remove("^\\.+") %>% str_remove("\\.+$") %>% str_replace_all("\\)", "\\(") %>% stringi::stri_reverse()) %>% 
  mutate(pistruct = case_when(((!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)"))) & wobble =="normal") ~ "perfect match",
                              ((!(str_detect(piRNAmatch,"\\.")) & !(str_detect(targetmatch,"\\)\\.+\\)"))) & wobble =="wobble") ~ "functional", # wobble pair
                              (str_detect(piRNAmatch,"^\\.\\(\\(\\(") & pilen ==(picount + 1)) &!(str_detect(targetmatch,"\\)\\.+\\)"))  ~ "functional", # firstbase mismatch
                              (str_detect(piRNAmatch,"\\(\\(\\(\\.$") & pilen ==(picount + 1)) &!(str_detect(targetmatch,"\\)\\.+\\)")) ~ "Miwi_tolerated", # lastbase mismatch
                              (piRNAmatch == target_rev & pilen ==(picount + 1))  ~ "Miwi_tolerated", # single mismatch
                              str_detect(piRNAmatch,"^\\.\\({19}") & str_detect(targetmatch, "\\){19}\\.+$") ~ "Miwi_tolerated", # 2-20bp is paired
                              str_detect(piRNAmatch,"^\\({20}") & str_detect(targetmatch, "\\){20}\\.+$") ~ "Miwi_tolerated", # 1-20bp is paired
                              TRUE ~ "not functional"))  %>% 
  filter(pistruct != "not functional")
Vienna1_model %>% write_csv(paste0(outdir, "piRNA_model.tsv"))
# ROC
Vienna2_all <- Vienna1_model %>% filter(RNAtype %in% c("TE","mRNA")) %>%  #filter(type3 !="rmst") %>%
  mutate(ki = exp((-dG * 4.184 * 1000)/(8.3145 * 299.15)),
         kiE  = ki * CPM) %>% 
  group_by(targetname, RNAtype,genesymbol) %>% 
  summarise(psi = sum(kiE),
            omega = sum (CPM),meanpilen = mean(pilen),
            meandG = mean(dG), pinum = n()) %>% ungroup() %>% 
  left_join(table3 %>% rename(targetname = name,Teyes = type) %>% select(targetname,ctrl_CPM,KD_CPM,logFC,Teyes) %>% 
              mutate(Teyes = Teyes %>% as.character()) %>% 
              mutate(Teyes = if_else(Teyes == "positive","independent",if_else(Teyes =="negative","independent",Teyes))), by = "targetname") %>% drop_na() %>% 
  mutate(phi = ctrl_CPM * psi, chi = phi / omega, log10CPM = log10(ctrl_CPM)) %>% 
  mutate(change = if_else(logFC>1, "change","stay"), tau = omega/ctrl_CPM) %>% 
  select(targetname, RNAtype, genesymbol,Teyes,psi,omega,phi,chi,tau,ctrl_CPM,KD_CPM, meandG,pinum,meanpilen,logFC)
Vienna2_all %>% write_csv(paste0(outdir, "piRNA_model2.csv"))
ROC_cal <- function(data, sorter,samplename) {
  sorter2 <- rlang::enquo(sorter)
  model_pi <- data %>% 
    arrange(desc(!!sorter2)) %>% 
    mutate(order = row_number() %>% as.integer(),
           tmp = if_else(Teyes == "target RNA",1,0),
           total_posi = sum(tmp),
           tmp2 = if_else(Teyes == "target RNA",0,1),
           total_nega = sum(tmp2),
           counter_posi = cumsum(tmp)/total_posi,
           counter_nega = cumsum(tmp2)/total_nega,
           sample = samplename) 
  return(model_pi)
}
model_psi <- ROC_cal(data = Vienna2_all ,sorter = psi,samplename = "psi")
model_omega <- ROC_cal(data = Vienna2_all ,sorter = omega,samplename = "omega")
model_phi <- ROC_cal(data = Vienna2_all ,sorter = phi,samplename = "phi")
model_chi <- ROC_cal(data = Vienna2_all ,sorter = chi,samplename = "chi")
model_tau <- ROC_cal(data = Vienna2_all ,sorter = tau,samplename = "tau")
model_meandG <- ROC_cal(data = Vienna2_all ,sorter = meandG,samplename = "meandG")
model_pinum <- ROC_cal(data = Vienna2_all ,sorter = pinum,samplename = "pinum")
model_meanpilen <- ROC_cal(data = Vienna2_all ,sorter = meanpilen,samplename = "meanpilen")

x_determines <- bind_rows(model_psi,model_omega,model_phi,model_chi,model_tau, model_meandG, model_pinum, model_meanpilen) %>%
  filter(counter_posi != 0 & counter_posi != 1) %>% 
  group_by(sample) %>% summarise(mean_x = mean(counter_nega), mean_y = mean(counter_posi)) %>% ungroup() 

bind_rows(model_psi,model_omega,model_phi,model_chi,model_tau,model_meandG, model_pinum, model_meanpilen) %>%
  group_by(sample) %>% mutate(mean_x = mean(counter_nega), mean_y = mean(counter_posi)) %>% ungroup() %>% 
  ggplot(aes(x = counter_nega, y = counter_posi, color = sample))+
  geom_line(alpha = 0.8, size=1)+ geom_abline(aes(intercept = 0, slope = 1), color = "gray40",linetype = "dashed") +
  geom_text(data = x_determines, aes(label = sample, x= mean_x, y=mean_y, color = sample), size = 2) + 
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1))+scale_x_continuous(limits = c(0, 1))+
  labs(x = "False negative ratio",y = "True positive ratio") +
  theme(legend.title = element_blank(),legend.position = c(0.8,0.4),panel.grid.minor = element_blank()) +
  ggsave(paste0(outdir,"model_ROC_2.png"), width =5, height =5, dpi = 300)+
  ggsave(paste0(outdir,"model_ROC_2.pdf"), width =5, height =5)

bind_rows(model_tau %>% mutate(samplename = "piRNA/mRNA"),
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
AUC_psi <- AUC_cal(data = model_psi)
AUC_omega <- AUC_cal(data = model_omega)
AUC_phi <- AUC_cal(data = model_phi)
AUC_chi <- AUC_cal(data = model_chi)
AUC_tau <- AUC_cal(data = model_tau)

AUC_table <- bind_rows(AUC_psi %>% mutate(samplename = "psi"),
                       AUC_phi %>% mutate(samplename = "phi"),
                       AUC_chi %>% mutate(samplename = "chi"),
                       AUC_omega%>% mutate(samplename = "piRNA expression"),
                       AUC_tau %>% mutate(samplename = "piRNA/mRNA rate")) %>% 
  write_csv(paste0(outdir,"AUC_table.csv"))

#######################################################################################################
# scatter plot
#######################################################################################################
#edit here#####################################################################################################
change_order <- c("target RNA","independent")
Vienna2_all %>% mutate(Teyes = Teyes %>% factor(levels = change_order)) %>% 
  ggplot(aes(x = omega, y =tau, color = Teyes))+
  geom_point(alpha = 0.8)+scale_x_log10()+scale_y_log10() + theme_minimal()+
  ggrepel::geom_text_repel(data = subset(Vienna2_all %>% filter(genesymbol %in% c("CG12717", "Ste:CG33237","Dmel-hopper2", "Dmel-412","Dmel-Circe","Dmel-transib2","Dmel-Doc"))),
                           aes(label =genesymbol), size = 2.7, color = "magenta4",  point.padding = 0.1,max.overlaps =Inf,
                           box.padding = unit(1.2, "lines"), min.segment.length = 0)+
  ggrepel::geom_text_repel(data = subset(Vienna2_all %>% filter(genesymbol %in% c("Dmel-Tabor","Dmel-mdg1","Dmel-ZAM","Dmel-flea"))),
                           aes(label =genesymbol), size = 2.7, color = "black",  point.padding = 0.1,max.overlaps =Inf,nudge_x = 0.5, nudge_y =-0.5,
                           box.padding = unit(1.2, "lines"), min.segment.length = 0)+
  labs(x = "functional piRNA's \ntotal expression (CPM)",y = "piRNA/mRNA (CPM)") +
  theme(legend.title = element_blank(),legend.position = "bottom")+
  ggsave(paste0(outdir,"TE_piCPM_rate2.png"), width =3.5, height =3.8, dpi = 300)+
  ggsave(paste0(outdir,"TE_piCPM_rate2.pdf"), width =3.5, height =3.8)

change_order2 <- c("target RNA","independent", "mRNA")
Vienna2_all %>% mutate(Teyes = if_else(Teyes=="independent" & RNAtype =="mRNA", "mRNA", Teyes) %>% factor(levels = change_order2)) %>% 
  ggplot(aes(x = omega, y =logFC, color = Teyes))+
  geom_point(alpha = 0.8)+ scale_x_log10()+ theme_minimal()+
  ggrepel::geom_text_repel(data = subset(Vienna2_all %>% filter(genesymbol %in% c("CG12717", "Ste:CG33237","Dmel-hopper2", "Dmel-412","Dmel-Circe","Dmel-transib2","Dmel-Doc"))),
                           aes(label =genesymbol), size = 2.7, color = "magenta4",  point.padding = 0.1,max.overlaps =Inf,
                           box.padding = unit(1.2, "lines"), min.segment.length = 0)+
  ggrepel::geom_text_repel(data = subset(Vienna2_all %>% filter(genesymbol %in% c("Dmel-Tabor","Dmel-mdg1","Dmel-ZAM","Dmel-flea"))),
                           aes(label =genesymbol), size = 2.7, color = "black",  point.padding = 0.1,max.overlaps =Inf,nudge_x = 0.5, nudge_y =-0.5,
                           box.padding = unit(1.2, "lines"), min.segment.length = 0)+
  labs(x = "piRNA expression (CPM)\n(perfect and near-perfect match)",y = "log2FC") +
  theme(legend.title = element_blank(),legend.position = "bottom")+
  ggsave(paste0(outdir,"TE_piCPM_WT2.png"), width =4, height =3.8, dpi = 300)+
  ggsave(paste0(outdir,"TE_piCPM_WT2.pdf"), width =4, height =3.8)
###############################################################################################################
