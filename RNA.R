
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gt)))
suppressMessages(suppressWarnings(require(edgeR)))
suppressMessages(suppressWarnings(require(ggrepel)))
suppressMessages(suppressWarnings(require(GGally)))
suppressMessages(suppressWarnings(require(ggthemes)))

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

outdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/"

#######################################################################################################
# CAGE count data detection
#######################################################################################################

CAGEdir <- "/media/hermione/CAGE/figures2/"
cage_anno2 <- read_tsv(paste0(CAGEdir,"ctss_annotation.tsv"))

file_cage <- dir(paste0(CAGEdir,"ctss/"), pattern = ".cov$", full.names = TRUE)
file_cage2 <- file_cage  %>%
  str_remove(".cov$") %>%
  str_remove(paste0(CAGEdir,"ctss/", "/")) %>%
  str_remove("cage")
cagedf <- vector("list", length(file_cage))
for (i in seq_along(file_cage)) {
  cagedf[[i]] <- read_tsv(file_cage[[i]],
                          col_names = c("chr1","start1","end1","strand1",
                                        "chr2","start2","end2","name","count","strand2"))
  cagedf[[i]][["sample"]] <- file_cage2[[i]]
}
cage_all <- cagedf %>% bind_rows() %>% filter(strand1 == strand2) %>%
  group_by(chr1, start1, end1, strand1, sample) %>% summarise(count = sum(count)) %>% ungroup()

cage2 <- cage_all %>% arrange(chr1, start1, end1) %>% unite(chr1, start1, end1,strand1, sep = "_", col = "name") %>% 
  left_join(cage_anno2 %>% unite(chr, start, end,direction, sep = "_", col = "name"), by = "name") %>% 
  drop_na()

# TE
CAGETEdir <- "/media/hermione/CAGE/TE/"
file_cage <- dir(CAGETEdir, pattern = "_TE.bed$", full.names = TRUE)
file_cage2 <- file_cage  %>%
  str_replace_all("_TE.bed$", "") %>%
  str_replace(paste0(CAGETEdir, "/"), "")
cagedf <- vector("list", length(file_cage))
for (i in seq_along(file_cage)) {
  cagedf[[i]] <- vroom::vroom(file_cage[[i]],
                              col_names = c("chr1","start1","end1","chr2","start2","end2","name","qual","strand1","strand2")) %>%
    mutate(point = if_else(strand1 == "+",start1+1, end1-1)) %>%
    group_by(chr1, point,strand1) %>% count()
  cagedf[[i]][["sample"]] <- file_cage2[[i]]
}
cage_TE <- cagedf %>% bind_rows() %>%
  pivot_wider(names_from = sample, values_from = n, values_fill =list(n=0)) %>% ungroup() %>%
  mutate(siEGFP_rep1 = siEGFPrep1_1+siEGFPrep1_2,
         siEGFP_rep2 = siEGFPrep2_1+siEGFPrep2_2,
         siEGFP_rep3 = siEGFPrep3_1+siEGFPrep3_2,
         siPIWI_rep1 = siPIWIrep1_1+siPIWIrep1_2,
         siPIWI_rep2 = siPIWIrep2_1+siPIWIrep2_2,
         siPIWI_rep3 = siPIWIrep3_1+siPIWIrep3_2) %>%
  select(chr1,point,strand1, siEGFP_rep1,siEGFP_rep2,siEGFP_rep3,siPIWI_rep1,siPIWI_rep2,siPIWI_rep3) %>%
  pivot_longer(-c(chr1,point, strand1), names_to = "sample", values_to = "count") %>%
  filter(count != 0) %>% filter(strand1 =="+") %>% rename(genename = chr1) %>% 
  mutate(name = paste0(genename,"_",point), genesymbol = genename, anno = genename, tfsymbol = genename) %>% 
  select(name, sample,count,genename,genesymbol, anno, tfsymbol)

cage_matrix <- bind_rows(cage2,cage_TE ) %>% 
  write_csv(paste0(outdir, "CAGE_Matrix.csv"))

# stats of CAGE mapping

cage_matrix <- read_csv(paste0(outdir, "CAGE_Matrix.csv"))

cage_matrix %>% filter(str_detect(genename, "FBgn")) %>% select(name) %>% distinct() %>% mutate(type = "total") %>% 
  group_by(type) %>% count() %>% write_csv(paste0(outdir, "CAGE_isoform_table1.csv"))
cage_matrix %>% filter(str_detect(genename, "FBgn")) %>% group_by(sample) %>% 
  summarise(isoform_num = n(), total_count = sum(count)) %>% write_csv(paste0(outdir, "CAGE_isoform_table2.csv"))

cage_matrix %>% filter(!(str_detect(genename, "FBgn"))) %>% group_by(sample) %>% 
  summarise(isoform_num = n(), total_count = sum(count)) %>% write_csv(paste0(outdir, "CAGE_TE_table1.csv"))
# creating normalization factor

cage_normalize_factor_table <- cage_matrix %>% 
  mutate(KD = sample %>% str_remove("_rep\\d+")) %>% 
  group_by(KD,sample) %>% summarise(total = sum(count)) %>% ungroup() %>% 
  group_by(KD) %>% mutate(totalKD = sum(total)) %>% ungroup()
  
cage_normalize_factor_table %>% 
  write_csv(paste0(outdir, "CAGE_nomalize_factor.csv"))

fusion_target <- c("Oat","CG42565","CG15209","CHKov1","CG2233","NT5E-2","CG8303","TwdlL","Adgf-A","Irk1","CG7460","Cyp12e1")
dependent_TE <- c("297","412","gypsy","gtwin","blood","17.6","mdg1","Tabor","I-element","HMS-Beagle2",
                  "Quasimodo", "Stalker","Stalker2","Stalker4","Idefix","mdg3","F-element","gypsy5","Juan")
change_order <- c("negative","positive","Piwi dependent TE","FDR top12")

cage_matrix_scatter <- cage_matrix %>% mutate(KD = sample %>% str_remove("_rep\\d+")) %>% 
  group_by(KD,genesymbol,genename) %>% summarise(count = sum(count)) %>% ungroup() %>% 
  group_by(KD) %>% mutate(CPM = 1000000 * count/ sum(count)) %>% ungroup() %>% select(-count) %>% 
  pivot_wider(names_from = "KD", values_from = "CPM", values_fill = 0) %>% 
  mutate(per =  log2((siPIWI+1)/(siEGFP+1)),
         type = if_else(abs(per) > 2 & (siEGFP + siPIWI) > 0.5, "positive", "negative")) %>% 
  mutate(type = case_when(genesymbol %in% paste0("Dmel-",dependent_TE) ~ "Piwi dependent TE",
                          genesymbol %in% fusion_target ~ "FDR top12",
                          TRUE ~ type)%>% factor(levels = change_order)) 

myfun <- function(x) log10(((exp(x*log(10))-1) -3)/4 +1)
myfun2 <- function(x) log10((4 * (exp(x*log(10))-1) + 3) +1)

cage_matrix_scatter %>%  
  ggplot(aes(x = log10(siEGFP + 1), y = log10(siPIWI + 1), color = type))+
  geom_function(fun = myfun, color = "gray60",linetype = "dashed", xlim=c(0.7,4))+
  geom_function(fun = myfun2, color = "gray60",linetype = "dashed", xlim=c(0,3.5))+
  geom_point(size = 1.4,alpha = 0.6)+ scale_color_manual(values = c( "#808080","#dc143c","blue","seagreen4")) +
  geom_text_repel(data = subset(cage_matrix_scatter, type == "FDR top12"),
                  aes(label = genesymbol), color = "darkgreen", size = 2.5,max.overlaps = Inf,
                  segment.color = "gray90", fontface = "bold", segment.size = 0.5,nudge_x = -0.1,nudge_y = 1.8)+
  geom_text_repel(data = subset(cage_matrix_scatter, type == "Piwi dependent TE"),
                  aes(label = genesymbol), color = "blue4", size = 2.5,max.overlaps = Inf,
                  segment.color = "gray90", fontface = "bold", segment.size = 0.5,nudge_x = -0.2,nudge_y = 0.05)+
# geom_text_repel(data = subset(cage_matrix_scatter, type == "positive"),
#                 aes(label = genesymbol), color = "black", size = 2.5,max.overlaps = Inf,
#                 segment.color = "gray40", fontface = "bold", segment.size = 0.5)+
  geom_text_repel(data = subset(cage_matrix_scatter, genesymbol == "piwi"),
                  aes(label = genesymbol), color = "black", size = 2.5,max.overlaps = Inf,
                  segment.color = "gray40", fontface = "bold", segment.size = 0.5)+
  labs(x = "log10(EGFP KD (CPM) + 1)",y = "log10(Piwi KD (CPM) + 1)")+
  theme_minimal()+ guides(color = "none") +
  ggsave(paste0(outdir,"CAGE_scatter.png"), width =5, height =5, dpi = 600)
# edgeR 
# .Machine$double.xmin <- 2.225074e-500
# https://support.bioconductor.org/p/95224/
# If edgeR turns p-value = 0, set this value.
cage_matrix %>% select(genename, sample,count) %>% 
  group_by(genename, sample) %>% summarise(count = sum(count)) %>% ungroup() %>% arrange(genename,sample) %>% 
  pivot_wider(names_from = "sample",values_from = "count" , values_fill =list(count=0)) %>% 
  write_csv(paste0(outdir, "CAGE_Matrix_gene.csv"))

cage_edgeR <- cage_matrix %>% select(genename, sample,count) %>% 
  group_by(genename, sample) %>% summarise(count = sum(count)) %>% ungroup() %>% arrange(genename,sample) %>% 
  pivot_wider(names_from = "sample",values_from = "count" , values_fill =list(count=0)) %>% 
  column_to_rownames(var = "genename") %>% 
  select(siEGFP_rep1, siEGFP_rep2,siEGFP_rep3,siPIWI_rep1,siPIWI_rep2,siPIWI_rep3) %>% 
  mutate(total = rowSums(across(everything())),
         total2 = sum(total),
         total_CPM = 1000000*total/total2) %>% 
  filter(total_CPM >= 0.5) %>% 
  select(-total, -total2,-total_CPM)

count <- cage_edgeR
group <- factor(c("siEGFP", "siEGFP","siEGFP","siPIWI","siPIWI","siPIWI"))
d <- DGEList(counts = count, group = group)
d2 <- calcNormFactors(d, method = "RLE")
design <- model.matrix( ~ 0 + group, data = d)
d3 <- estimateGLMCommonDisp(d2, design)
result <- exactTest(d3)
table <- as.data.frame(topTags(result, n = nrow(count))) %>%
  rownames_to_column(var = "name") %>% as_tibble() %>%
  mutate(fdr2 = ifelse(FDR < 0.0001, "change", "stay"))

table2 <- table %>% left_join(cage_matrix %>% select(genename,genesymbol) %>% distinct() %>% rename(name = genename), by = "name") %>%
  mutate(logp = -log10(PValue), logFDR = -log10(FDR),
         state = case_when(logp > 50 & logFC > 1 ~  "plus change",
                           logp > 50 & logFC < -1 ~  "minus change",
                           TRUE ~ "stay"))

ggplot(table, mapping = aes(x = logCPM, y= logFC, colour = fdr2)) +
  geom_point(size = 1,alpha = 0.5) +
  scale_color_manual(values = c("#dc143c", "#808080")) +
  ggrepel::geom_text_repel(data = subset(table2 %>% filter(state %in% c("plus change", "minus change"))),
                           aes(label =genesymbol), size = 2.5, segment.color = "gray26",fontface = "bold",
                           segment.size = 0.5, color = "black",max.overlaps =Inf) +
  xlab("average logCPM") + ylab("log2FC") +
  guides(color = "none")+
  theme_minimal() +
  ggsave(paste0(outdir,"CAGE_edgeR_RLE.png"), width =8, height =8, dpi = 600)

table2 %>% 
  write_csv(paste0(outdir, "CAGE_edgeR.csv"))

png(paste0(outdir,"CAGE_corr.png"), width = 4*ppi, height = 4*ppi, res = ppi)
ggpairs(data = cage_edgeR %>% rename(`Piwi KD rep1` = siPIWI_rep1,`Piwi KD rep2` = siPIWI_rep2,`Piwi KD rep3` = siPIWI_rep3,
                                     `EGFP KD rep1` = siEGFP_rep1,`EGFP KD rep2` = siEGFP_rep2,`EGFP KD rep3` = siEGFP_rep3) %>% 
          mutate_at(vars(contains("rep")),list(~log10((. * 1000000/sum(.))+1))),
              xlab = "log10(CPM + 1) ",ylab = "log10(CPM + 1)",
              lower = list(continuous = wrap("points", alpha = 0.2, size = 0.1)),
              upper = list(continuous = wrap(cor_func, method = "pearson", symbol = expression('p ='))))+
  theme(text = element_text(size = 6), panel.grid.minor = element_blank())
dev.off()

# dm6 TSS change
cage_matrix <- read_csv(paste0(outdir, "CAGE_Matrix.csv")) %>% filter(str_detect(genename, "FBgn")) %>% 
  mutate(KD = sample %>% str_remove("_rep\\d+")) %>% 
  group_by(KD,genesymbol,genename,name) %>% summarise(count = sum(count)) %>% ungroup() %>% 
  group_by(KD) %>% mutate(CPM = 1000000 * count/ sum(count)) %>% ungroup() %>% select(-count) %>% 
  pivot_wider(names_from = "KD", values_from = "CPM", values_fill = 0)
table2 <- read_csv(paste0(outdir, "CAGE_edgeR.csv")) %>% 
  mutate(fdr2 = case_when(fdr2 == "change" & logFC > 0 ~ "up",
                          fdr2 == "change" & logFC < 0 ~ "down",
                          TRUE ~ fdr2))

TSS_change <- cage_matrix %>% 
  left_join(table2 %>% mutate(genename = name) %>% select(genename, fdr2), by = "genename") %>% drop_na() %>% 
  filter(siEGFP + siPIWI > 1) %>% filter(fdr2 != "stay") %>% 
  filter((siPIWI/siEGFP > 2 & fdr2 == "up") | (siPIWI/siEGFP < 0.5 & fdr2 == "down")) %>% 
  separate(name, sep = "_", into = c("chr","start","end","strand"))

TSS_change %>% mutate(name = paste0(fdr2, "@",genesymbol), dif = log2(siPIWI/siEGFP)) %>% filter(fdr2 == "up") %>% 
  select(chr,start,end,name, dif,strand) %>% write_tsv(paste0(outdir, "CAGE_genome_change_up.bed"), col_names = FALSE)
TSS_change %>% mutate(name = paste0(fdr2, "@",genesymbol), dif = log2(siPIWI/siEGFP)) %>% filter(fdr2 == "down") %>% 
  select(chr,start,end,name, dif,strand) %>% write_tsv(paste0(outdir, "CAGE_genome_change_down.bed"), col_names = FALSE)


#######################################################################################################
# Iso-seq count data detection
#######################################################################################################
RNAtype <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_RNAtype.tsv")
RNAtype2 <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_use.tsv") %>% 
  filter(!(type %in% c("gene","CDS","exon","start_codon","stop_codon","5UTR","3UTR"))) %>% 
  select(name, genesymbol,transcript)

cov_gene <- read_tsv("/media/ariel/Isoseq/merge/gene.collapsed.filtered_classification.txt")
# Extended Data, isoform length

cov_gene %>% mutate(fl = FL.siEGFP + FL.siPIWI,fl2 =  FL.siEGFP2 + FL.siPIWI2) %>% 
  mutate(length = 100*(length%/%100) ) %>% mutate(length = length/1000) %>%  group_by(length) %>%
  summarise(rep1 = sum(fl),rep2 = sum(fl2)) %>% ungroup() %>% 
  pivot_longer(c(rep1,rep2), names_to = "dataset",values_to = "count") %>% filter(count != 0) %>% 
  ggplot(aes(x = length, y = count,fill = dataset))+
  geom_col(alpha = 0.5 ,position = 'identity')+ 
  scale_y_log10() + scale_x_continuous(breaks = c(0,1,2,3,4,5,7.5,10)) +
  theme_minimal()+ theme(panel.grid.minor= element_blank(), legend.title = element_blank(),legend.position = c(0.8,0.8)) +
  labs(x = "transcript length (kbp)", y = "total fl count") +
  ggsave(paste0(outdir,"Isoseq_transcript_len_qc2_dataset.png"), width =4, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"Isoseq_transcript_len_qc2_dataset.pdf"), width =4, height =2.5)


cov_gene2 <- cov_gene %>% 
  filter(structural_category %in%  c("full-splice_match","incomplete-splice_match","novel_in_catalog", "novel_not_in_catalog")) %>% 
  left_join(RNAtype %>% select(name, genesymbol) %>% rename(associated_gene = name), by = "associated_gene" ) %>% 
  left_join(RNAtype2 %>% rename(associated_transcript = transcript, genesymbol2 = genesymbol), by = "associated_transcript" ) %>% 
  mutate(associated_gene = if_else(is.na(name), associated_gene, name),
         genesymbol = if_else(is.na(genesymbol),genesymbol2, genesymbol)) %>% 
  select(-name,-genesymbol2) %>% 
  group_by(associated_gene) %>% mutate(fl_EGFP = sum(FL.siEGFP+FL.siEGFP2),fl_PIWI = sum(FL.siPIWI+FL.siPIWI2)) %>% 
  ungroup() %>% mutate(per = log2((fl_PIWI+1)/(fl_EGFP+1))) %>% 
  select(-n_indels,-n_indels_junc,-FL,-dist_to_polyA_site, -within_polyA_site) 

# TE count
cov_TE <- read_tsv("/media/ariel/Isoseq/mergeTE/TE.collapsed.filtered_classification.txt")
cov_TE2 <- cov_TE %>% 
  filter(strand == "+") %>% 
  group_by(chrom) %>% mutate(fl_EGFP = sum(FL.siEGFP+FL.siEGFP2),fl_PIWI = sum(FL.siPIWI+FL.siPIWI2)) %>% 
  ungroup() %>% mutate(per = log2((fl_PIWI+1)/(fl_EGFP+1))) %>% 
  select(-n_indels,-n_indels_junc,-FL,-dist_to_polyA_site, -within_polyA_site) 

cov_Isoseq <- bind_rows(cov_gene2 %>% select(associated_gene, genesymbol, fl_EGFP,fl_PIWI,per) %>% distinct(),
                        cov_TE2 %>% mutate(associated_gene = chrom, genesymbol = chrom) %>% 
                          select(associated_gene, genesymbol, fl_EGFP,fl_PIWI,per) %>% distinct())

cov_Isoseq %>% 
  write_csv(paste0(outdir, "Isoseq_Matrix.csv"))


# stats of Isoseq mapping

cov_gene <- read_tsv("/media/ariel/Isoseq/merge/gene.collapsed.filtered_classification.txt")
cov_gene %>% group_by(structural_category) %>% 
  summarise(count = n(),FL = sum(FL.siEGFP,na.rm = TRUE) +sum(FL.siEGFP2,na.rm = TRUE) +
              sum(FL.siPIWI,na.rm = TRUE) + sum(FL.siPIWI2,na.rm = TRUE),
            FL_EGFP = sum(FL.siEGFP,na.rm = TRUE)  + sum(FL.siEGFP2,na.rm = TRUE),
            FL_PIWI = sum(FL.siPIWI,na.rm = TRUE)+ sum(FL.siPIWI2,na.rm = TRUE),
            FL_EGFP_rep1 = sum(FL.siEGFP,na.rm = TRUE), FL_EGFP_rep2 = sum(FL.siEGFP2,na.rm = TRUE),
            FL_PIWIrep1 = sum(FL.siPIWI,na.rm = TRUE), FL_PIWIrep2 =sum(FL.siPIWI2,na.rm = TRUE)) %>%
  write_csv(paste0(outdir, "Isoseq_isoform_table1.csv"))
cov_TE <- read_tsv("/media/ariel/Isoseq/mergeTE/TE.collapsed.filtered_classification.txt")
cov_TE %>% group_by(structural_category) %>% 
  summarise(count = n(),FL = sum(FL.siEGFP,na.rm = TRUE) +sum(FL.siEGFP2,na.rm = TRUE) +
              sum(FL.siPIWI,na.rm = TRUE) + sum(FL.siPIWI2,na.rm = TRUE),
            FL_EGFP = sum(FL.siEGFP,na.rm = TRUE) + sum(FL.siEGFP2,na.rm = TRUE),
            FL_PIWI = sum(FL.siPIWI,na.rm = TRUE)+ sum(FL.siPIWI2,na.rm = TRUE),
            FL_EGFP_rep1 = sum(FL.siEGFP,na.rm = TRUE), FL_EGFP_rep2 = sum(FL.siEGFP2,na.rm = TRUE),
            FL_PIWIrep1 = sum(FL.siPIWI,na.rm = TRUE), FL_PIWIrep2 =sum(FL.siPIWI2,na.rm = TRUE)) %>%
  write_csv(paste0(outdir, "Isoseq_TE_table1.csv"))

# creating normalization factor
Isoseq_normalize_factor_table <- cov_Isoseq %>% 
  pivot_longer(cols = c(fl_EGFP, fl_PIWI), names_to = "KD",values_to = "count") %>% 
  group_by(KD) %>% summarise(total = sum(count)) %>% ungroup() 

Isoseq_normalize_factor_table %>% 
  write_csv(paste0(outdir, "Isoseq_nomalize_factor.csv"))

# cor plot
cov_Isoseq_cor <- bind_rows( cov_gene %>% 
                           filter(structural_category %in%  c("full-splice_match","incomplete-splice_match","novel_in_catalog", "novel_not_in_catalog")) %>% 
                           left_join(RNAtype %>% select(name, genesymbol) %>% rename(associated_gene = name), by = "associated_gene" ) %>% 
                           left_join(RNAtype2 %>% rename(associated_transcript = transcript, genesymbol2 = genesymbol), by = "associated_transcript" ) %>% 
                           mutate(associated_gene = if_else(is.na(name), associated_gene, name),
                                  genesymbol = if_else(is.na(genesymbol),genesymbol2, genesymbol)) %>% 
                            rename(fl_EGFP1 = FL.siEGFP,fl_EGFP2 = FL.siEGFP2,fl_PIWI1 = FL.siPIWI,fl_PIWI2 = FL.siPIWI2) %>% 
                             select(associated_gene, genesymbol, fl_EGFP1,fl_PIWI1,fl_EGFP2,fl_PIWI2) %>% distinct(),
    cov_TE %>%  filter(strand == "+") %>%
     rename(fl_EGFP1 = FL.siEGFP,fl_EGFP2 = FL.siEGFP2,fl_PIWI1 = FL.siPIWI,fl_PIWI2 = FL.siPIWI2,associated_gene = chrom, genesymbol = chrom) %>% 
                          select(associated_gene, genesymbol, fl_EGFP1,fl_PIWI1,fl_EGFP2,fl_PIWI2) %>% distinct()) %>% 
  group_by(associated_gene, genesymbol) %>% summarise(siEGFP_rep1 = sum(fl_EGFP1),siEGFP_rep2 = sum(fl_EGFP2),
                                                      siPIWI_rep1 = sum(fl_PIWI1), siPIWI_rep2 = sum(fl_PIWI2)) %>% ungroup()

png(paste0(outdir,"Isoseq_corr.png"), width = 3*ppi, height = 3*ppi, res = ppi)
ggpairs(data = cov_Isoseq_cor %>%
          rename(`Piwi KD rep1` = siPIWI_rep1,`Piwi KD rep2` = siPIWI_rep2,`EGFP KD rep1` = siEGFP_rep1,`EGFP KD rep2` = siEGFP_rep2) %>% 
          mutate_at(vars(contains("rep")),list(~log10((. * 1000000/sum(.))+1))) %>% select(-associated_gene, -genesymbol),
        xlab = "log10(CPM + 1) ",ylab = "log10(CPM + 1)",
        lower = list(continuous = wrap("points", alpha = 0.2, size = 0.1)),
        upper = list(continuous = wrap(cor_func, method = "pearson", symbol = expression('p ='))))+
  theme(text = element_text(size = 6), panel.grid.minor = element_blank())
dev.off()
#

cov_Isoseq2 <- cov_Isoseq %>% 
  mutate(EGFP_CPM = 1000000 * fl_EGFP/ sum(fl_EGFP),
         PIWI_CPM = 1000000 * fl_PIWI/ sum(fl_PIWI),
         per =  log2((PIWI_CPM+1)/(EGFP_CPM+1)),
         type = if_else(abs(per) > 2 & (fl_EGFP + fl_PIWI) > 20, "positive", "negative"))  %>% 
  mutate(type = case_when(genesymbol %in% paste0("Dmel-",dependent_TE) ~ "Piwi dependent TE",
                          genesymbol %in% fusion_target ~ "FDR top12",
                          TRUE ~ type)%>% factor(levels = change_order)) 
myfun <- function(x) log10(((exp(x*log(10))-1) -3)/4 +1)
myfun2 <- function(x) log10((4 * (exp(x*log(10))-1) + 3) +1)

cov_Isoseq2 %>%  
  ggplot(aes(x = log10(EGFP_CPM + 1), y = log10(PIWI_CPM + 1), color = type))+
  geom_function(fun = myfun, color = "gray60",linetype = "dashed", xlim=c(0.7,4))+
  geom_function(fun = myfun2, color = "gray60",linetype = "dashed", xlim=c(0,3.5))+
  geom_point(size = 1.4,alpha = 0.6)+ scale_color_manual(values = c( "#808080","#dc143c","blue","seagreen4")) +
  geom_text_repel(data = subset(cov_Isoseq2, type == "FDR top12"),
                  aes(label = genesymbol), color = "darkgreen", size = 2.5,max.overlaps = Inf,
                  segment.color = "gray90", fontface = "bold", segment.size = 0.5,nudge_x = -0.15,nudge_y = 1.8)+
  geom_text_repel(data = subset(cov_Isoseq2, type == "Piwi dependent TE"),
                  aes(label = genesymbol), color = "blue4", size = 2.5,max.overlaps = Inf,
                  segment.color = "gray90", fontface = "bold", segment.size = 0.5,nudge_x = 0.08,nudge_y = -0.55)+
#  geom_text_repel(data = subset(cov_Isoseq2, type == "positive"),
#                  aes(label = genesymbol), color = "black", size = 2.5,max.overlaps = Inf,
#                  segment.color = "gray40", fontface = "bold", segment.size = 0.5)+
  geom_text_repel(data = subset(cov_Isoseq2, genesymbol == "piwi"),
                  aes(label = genesymbol), color = "black", size = 2.5,max.overlaps = Inf,
                  segment.color = "gray40", fontface = "bold", segment.size = 0.5)+
  labs(x = "log10(EGFP KD (CPM) + 1)",y = "log10(Piwi KD (CPM) + 1)")+
  theme_minimal()+ guides(color = "none") +
  ggsave(paste0(outdir,"Isoseq_scatter.png"), width =5, height =5, dpi = 600)



#######################################################################################################
# CAGE fusion count
#######################################################################################################

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
  summarise(count = sum(count)) %>% ungroup() 

chimeraall2 <- chimeraall %>% 
  pivot_wider(names_from = "sample",values_from = "count", values_fill =list(count=0))

chimera_write <- chimeraall %>% 
  filter(fusionclass == "valid") %>%
  write_csv(paste0(outdir,"CAGE_all_chimera_list.csv"))

# stats of CAGE chimera

chimera_write <- read_csv(paste0(outdir,"CAGE_all_chimera_list.csv"))

chimera_write %>% group_by(sample) %>% summarise(isoform = n(), count = sum(count)) %>% 
  write_csv(paste0(outdir, "CAGE_chimera_table1.csv"))
chimera_write %>% mutate(mRNAclass = if_else(str_detect(chr_donorA, "^chr"), "genome-TE","TE-genome")) %>% 
  group_by(sample,mRNAclass) %>% summarise(isoform = n(), count = sum(count)) %>% 
  write_csv(paste0(outdir, "CAGE_chimera_table2.csv"))

#######################################################################################################
# Iso-seq fusion count
#######################################################################################################
fl_count <- read_tsv("/media/ariel/Isoseq/fusion/polish_isoforms_fusion.fl_count.tsv")
merge_chimera <- read_tsv("/media/ariel/Isoseq/fusion/chimera/polish_isoforms_fusion_classification.txt") %>% 
  select(isoform,chrom,strand,length, exons,structural_category,associated_gene,associated_transcript,
         subcategory,RTS_stage,FSM_class,iso_exp,gene_exp,ratio_exp) %>%
  mutate(pbid = isoform %>% str_extract("PBfusion.\\d+"),
         isonum = isoform %>% str_extract("\\.\\d+$") %>% str_remove("\\.") %>% as.integer(),
         index = paste0("index",row_number()),
         associated_gene =  if_else(str_detect(associated_gene,"novelGene_\\d+"),
                                    str_replace_all(associated_gene,"novelGene_","novelGene"),associated_gene) %>% str_replace("_AS","@AS")) %>% 
  separate_rows(associated_gene,sep = "_") %>%
  left_join(RNAtype %>% rename(associated_gene = name), by = "associated_gene") %>% 
  mutate(type = if_else(is.na(type) & str_detect(chrom,"chr"), "novel", if_else(is.na(type), "TE",type)),
         genesymbol = if_else(is.na(genesymbol) & str_detect(chrom,"chr"), associated_gene, if_else(is.na(genesymbol), chrom,genesymbol))) %>%
  group_by(index,isoform,chrom,strand,length, exons,structural_category,associated_transcript,
           subcategory,RTS_stage,FSM_class,iso_exp,gene_exp,ratio_exp,pbid,isonum) %>%
  summarise(associated_gene = paste0(associated_gene, collapse = "_"),
            RNAtype = paste0(type, collapse = "_"),
            symbol = paste0(genesymbol, collapse = "_")) %>% ungroup() %>%
  left_join(fl_count %>% rename(pbid = id), by = "pbid") %>% 
  mutate(associated_gene = if_else(chrom %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY"), associated_gene,chrom),
         type = if_else(chrom %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY"), "genome","TE")) %>% 
  select(associated_gene,isoform,structural_category,associated_transcript,type,pbid,isonum,
         iso_exp,gene_exp,ratio_exp,RNAtype,symbol,siEGFP,siPIWI,siEGFP2,siPIWI2) %>%
  arrange(pbid,isonum) %>% group_by(pbid) %>% mutate(state = paste0(type,collapse = "-")) %>% ungroup() %>% 
  group_by(pbid,state,siEGFP,siPIWI,siEGFP2,siPIWI2) %>%
  summarise(first = associated_gene[1], firstsymbol = symbol[1],
            second = paste0(associated_gene %>% tail(-1), collapse = "@"),
            secondsymbol = paste0(symbol %>% tail(-1), collapse = "@"),firstRNAtype = RNAtype[1],
            secondRNAtype = paste0(RNAtype %>% tail(-1), collapse = "@")) %>% ungroup() %>% 
  mutate(EGFP = siEGFP + siEGFP2, PIWI = siPIWI + siPIWI2) %>% select(-c(siEGFP, siPIWI,siEGFP2, siPIWI2)) %>% 
  pivot_longer(cols = c("EGFP","PIWI"), values_to = "count_fl", names_to = "KD") 
  
merge_chimera2 <- merge_chimera %>% group_by(KD,first,second,state,firstsymbol,secondsymbol,firstRNAtype,secondRNAtype) %>%
  summarise(pbid = paste0(pbid,collapse = "@"), total_fl = sum(count_fl), total_iso = n()) %>% ungroup() %>%
  mutate(chimera = paste0(first,"&",second),
         chimerasymbol = paste0(firstsymbol,"&",secondsymbol),
         chimeraRNAtype = paste0(firstRNAtype,"&",secondRNAtype)) 

merge_chimera %>% write_csv(paste0(outdir,"Isoseq_chimera_list.csv"))

  
# stats of Isoseq mapping
merge_chimera_table <- read_tsv("/media/ariel/Isoseq/fusion/chimera/polish_isoforms_fusion_classification.txt") %>% 
  select(isoform,chrom,strand,length, exons,structural_category,associated_gene,associated_transcript,
         subcategory,RTS_stage,FSM_class,iso_exp,gene_exp,ratio_exp) %>%
  mutate(pbid = isoform %>% str_extract("PBfusion.\\d+"),
         isonum = isoform %>% str_extract("\\.\\d+$") %>% str_remove("\\.") %>% as.integer(),
         index = paste0("index",row_number()),
         associated_gene =  if_else(str_detect(associated_gene,"novelGene_\\d+"),
                                    str_replace_all(associated_gene,"novelGene_","novelGene"),associated_gene) %>% str_replace("_AS","@AS")) %>% 
  separate_rows(associated_gene,sep = "_") %>%
  left_join(RNAtype %>% rename(associated_gene = name), by = "associated_gene") %>% 
  mutate(type = if_else(is.na(type) & str_detect(chrom,"chr"), "novel", if_else(is.na(type), "TE",type)),
         genesymbol = if_else(is.na(genesymbol) & str_detect(chrom,"chr"), associated_gene, if_else(is.na(genesymbol), chrom,genesymbol))) %>%
  group_by(index,isoform,chrom,strand,length, exons,structural_category,associated_transcript,
           subcategory,RTS_stage,FSM_class,iso_exp,gene_exp,ratio_exp,pbid,isonum) %>%
  summarise(associated_gene = paste0(associated_gene, collapse = "_"),
            RNAtype = paste0(type, collapse = "_"),
            symbol = paste0(genesymbol, collapse = "_")) %>% ungroup() %>%
  left_join(fl_count %>% rename(pbid = id), by = "pbid") %>% 
  mutate(associated_gene = if_else(chrom %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY"), associated_gene,chrom),
         type = if_else(chrom %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY"), "genome","TE")) %>% 
  select(associated_gene,isoform,structural_category,associated_transcript,type,pbid,isonum,
         iso_exp,gene_exp,ratio_exp,RNAtype,symbol,siEGFP,siPIWI,siEGFP2,siPIWI2) %>%
  arrange(pbid,isonum) %>% group_by(pbid) %>% mutate(state = paste0(type,collapse = "-")) %>% ungroup() %>% 
  group_by(pbid,state,siEGFP,siPIWI,siEGFP2,siPIWI2) %>%
  summarise(first = associated_gene[1], firstsymbol = symbol[1],
            second = paste0(associated_gene %>% tail(-1), collapse = "@"),
            secondsymbol = paste0(symbol %>% tail(-1), collapse = "@"),firstRNAtype = RNAtype[1],
            secondRNAtype = paste0(RNAtype %>% tail(-1), collapse = "@")) %>% ungroup()

merge_chimera_table %>% filter(str_detect(state, "TE")) %>% filter(str_detect(state, "genome")) %>% select(pbid) %>% distinct() %>% mutate(type = "total") %>% 
  group_by(type) %>% count() %>% write_csv(paste0(outdir, "Isoseq_chimera_table1.csv"))
merge_chimera_table %>% pivot_longer(cols = c("siEGFP","siPIWI","siEGFP2","siPIWI2"), values_to = "count_fl", names_to = "KD") %>% 
  filter(str_detect(state, "TE")) %>% filter(str_detect(state, "genome")) %>% filter(count_fl!=0) %>% group_by(KD) %>% 
  summarise(isoform_num = n(), total_count = sum(count_fl)) %>% write_csv(paste0(outdir, "Isoseq_chimera_table2.csv"))

# creating normalization factor
merge_chimera <- read_csv(paste0(outdir,"Isoseq_chimera_list.csv"))
dependent_TE <- c("297","412","gypsy","gtwin","blood","17.6","mdg1","Tabor","I-element","HMS-Beagle2",
                  "Quasimodo", "Stalker","Stalker2","Stalker4","Idefix","mdg3","F-element","gypsy5","Juan")

merge_chimera_scatter <- merge_chimera2  %>% filter(str_detect(state, "TE")) %>% filter(str_detect(state, "genome")) %>% 
  select(chimera, KD,total_fl, chimerasymbol,chimeraRNAtype) %>%
  pivot_wider(names_from = KD,values_from = total_fl, values_fill =list(total_fl=0)) %>% 
  mutate(relatedTE = chimera %>%
           str_extract_all("Dmel-17.6&|Dmel-17.6@|Dmel-17.6$|Dmel-\\w+&|Dmel-\\w+@|Dmel-\\w+$|Dmel-\\w+-\\w+&|Dmel-\\w+-\\w+@|Dmel-\\w+-\\w+$") %>% 
           str_remove_all("&|@")) 
  
merge_chimera_scatter2 <- merge_chimera_scatter %>% 
  filter(str_detect(relatedTE,"c\\(")) %>% 
  separate(relatedTE, sep = " ", into = c("relatedTE1","relatedTE2")) %>% 
  mutate(relatedTE1 = relatedTE1 %>% str_remove("c\\(") %>% str_remove_all('"') %>% str_remove(","),
         relatedTE2 = relatedTE2 %>% str_remove("\\)") %>% str_remove_all('"'))
merge_chimera_scatter3 <- merge_chimera_scatter %>% 
  left_join(merge_chimera_scatter2 %>% select(chimera,EGFP,PIWI,relatedTE1,relatedTE2), by = c("chimera","EGFP","PIWI")) %>%
  filter(!(str_detect(chimera,"novelGene"))) %>% 
  mutate(type = case_when(relatedTE %in% paste0("Dmel-",dependent_TE) ~ "Piwi dependent TE",
                          relatedTE1 %in% paste0("Dmel-",dependent_TE) ~ "Piwi dependent TE",
                          relatedTE2 %in% paste0("Dmel-",dependent_TE) ~ "Piwi dependent TE",
                          TRUE ~ "Piwi independent TE")) %>% 
  mutate(FC = log2((PIWI+1)/(EGFP +1)))

# percentage of TEs
fusion_overlap <- function(overlapdir, name, data) {
  overlap <- read_tsv(overlapdir,
                      col_names = c("chr1","start1","end1","name1","score1","strand1","chr2","source","type",
                                    "gffstart","gffend","score2","strand2","score3","name2"))
  overlap2 <- overlap %>%
    separate(name2,sep = ";",into = c("pbid","geneid","fil")) %>% 
    mutate(geneid = geneid %>% str_remove('gene_id "') %>% str_remove('"'),
           pbid = pbid %>% str_remove(' transcript_id "') %>% str_remove('"') %>% str_extract("PBfusion.\\d+")) %>%
    filter(type == "transcript")
  merge <- merge_chimera %>% filter(str_detect(firstRNAtype,"TE") | str_detect(secondRNAtype,"TE")) %>% filter(state %in% c("TE-genome", "genome-TE")) %>%
    mutate(relateTE = if_else(firstRNAtype=="TE", first, second) %>% str_remove_all("FBgn\\d+") %>% str_remove_all("novelGene\\d+") %>% 
             str_remove_all("novelGene_") %>% str_remove_all("@")) %>% 
    left_join(overlap2 %>% rename(RNAtype2 = type), by = "pbid") %>% 
    mutate(queryTE = name1) %>% separate(queryTE, sep = "@",into = c("queryTE","caller","fullness")) %>%
    filter(queryTE == relateTE)
  merge2 <- merge %>%
    mutate(distance1 = abs(gffstart-start1),
           distance2 = abs(gffend -start1),
           nearst = if_else(distance1 > distance2, distance2, distance1)) %>%
    select(pbid, relateTE, chr1,start1,end1,strand1,name1,nearst) %>% 
    group_by(pbid,relateTE) %>% filter(nearst == min(nearst)) %>% ungroup() %>%
    left_join(merge, by = c("pbid", "relateTE", "chr1","start1","end1","strand1","name1")) %>% distinct()
  return(merge2)
}

merge2_chimera <-fusion_overlap(overlapdir = "/media/pericles/TEfind/overlap/polish_overlap.txt",data = merge_chimera) %>%
  filter(str_detect(state, "TE")) %>% filter(str_detect(state, "genome"))
merge2_chimera %>%
  write_csv(paste0(outdir,"Isoseq_chimera_list_overlap.csv"))

#######################################################################################################
# CAGE and Iso-seq merge
#######################################################################################################
dependent_TE <- c("297","412","gypsy","gtwin","blood","17.6","mdg1","Tabor","I-element","HMS-Beagle2",
                  "Quasimodo", "Stalker","Stalker2","Stalker4","Idefix","mdg3","F-element","gypsy5","Juan")

bind_rows(merge_chimera_scatter3 %>% rename(siEGFP = EGFP,siPIWI=PIWI) %>% select(type, siEGFP,siPIWI) %>% mutate(seqmethod = "Iso-seq"),
          CAGE_chimera_scatter %>% select(type, siEGFP,siPIWI) %>% mutate(seqmethod = "CAGE")) %>% 
  filter((siEGFP + siPIWI) > 2) %>% 
  mutate(type = if_else(type == "Piwi dependent TE", "Piwi dependent TE-mRNA chimera","Piwi independent TE-mRNA chimera")) %>% 
  ggplot(aes(x = log2(siEGFP+1),y = log2(siPIWI+1)))+
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(alpha = 0.5, aes(color = type))+
  geom_density2d( color = "gray")+
  facet_grid(seqmethod~type) +
  theme_minimal()+
  guides(color ="none")+ labs(x = "log2(EGFP KD (count) + 1)",y = "log2(Piwi KD (count) + 1)")+
  ggsave(paste0(outdir,"merged_fusion_scatter.png"), width =5, height =5, dpi = 600)+
  ggsave(paste0(outdir,"merged_fusion_scatter.pdf"), width =5, height =5)

CAGE_TE_merge <- CAGE_chimera_scatter %>%
  mutate(TE = TEtype) %>%
  pivot_longer(cols = c("siEGFP", "siPIWI"), names_to = "KD", values_to = "count") %>% 
  left_join(cage_normalize_factor_table %>% select(KD,totalKD) %>% distinct(), by = "KD") %>% 
  mutate(CPM = 1000000 * count/totalKD) %>% 
  group_by(KD,TE) %>%
  summarise(CPM = sum(CPM)) %>% ungroup()
merge2_chimera <- read_csv(paste0(outdir,"Isoseq_chimera_list_overlap.csv"))
Isoseq_TE_merge <- merge2_chimera %>%
  mutate(TE = relateTE) %>%
  left_join(Isoseq_normalize_factor_table %>% mutate(KD = KD %>% str_remove("fl_")), by = "KD") %>% 
  mutate(CPM = 1000000 * count_fl/total) %>% 
  group_by(KD,TE) %>%
  summarise(CPM = sum(CPM)) %>% ungroup()

All_TE_merge <- CAGE_TE_merge %>% 
  pivot_wider(names_from = "KD", values_from = "CPM", values_fill = 0) %>% 
  rename(CAGE_siEGFP=siEGFP, CAGE_siPIWI=siPIWI) %>% 
  full_join(Isoseq_TE_merge %>% pivot_wider(names_from = "KD", values_from = "CPM", values_fill = 0) %>% 
              rename(Isoseq_siEGFP=EGFP, Isoseq_siPIWI=PIWI), by = "TE")
All_TE_merge %>% write_csv(paste0(outdir,"chimera_comp.csv"))

merge_chimera_order2 <- c("Dmel-gypsy","Dmel-gypsy5","Dmel-blood","Dmel-Juan","Dmel-17.6","Dmel-297","Dmel-412","Dmel-mdg1","Dmel-mdg3",
                          "Dmel-Stalker", "Dmel-Stalker2","Dmel-Stalker4","Dmel-Quasimodo","Dmel-gtwin","Dmel-I-element","Dmel-F-element","Dmel-Tabor", 
                         "Dmel-Idefix", "Dmel-HMS-Beagle2", "Dmel-springer","Dmel-copia","Dmel-roo","Dmel-jockey","Dmel-flea","Dmel-gypsy3","Dmel-Circe","Dmel-invader6", "others")

All_TE_merge2 <- All_TE_merge %>%
  mutate(TE = if_else(TE %in% merge_chimera_order2, TE,"others")) %>%
  pivot_longer(c(CAGE_siEGFP,CAGE_siPIWI,Isoseq_siEGFP,Isoseq_siPIWI), names_to = "KD",values_to = "CPM") %>% drop_na() %>% 
  group_by(KD,TE) %>%
  summarise(CPM = sum(CPM)) %>% ungroup() %>% 
  pivot_wider(names_from = "KD", values_from = "CPM", values_fill = 0) %>% 
  pivot_longer(c(CAGE_siEGFP,CAGE_siPIWI,Isoseq_siEGFP,Isoseq_siPIWI), names_to = "KD",values_to = "CPM") %>% 
  mutate(type = case_when(TE %in% paste0("Dmel-",dependent_TE) ~ "dependent",
                          TRUE ~ "independent"),
         my_col = if_else(type == "dependent", "red","black"),
         TE = TE %>% factor(levels = merge_chimera_order2),
         KD = KD %>% paste0(.," KD") %>% str_replace("PIWI","Piwi") %>% str_replace("_"," ") %>% 
           str_remove("si") %>% str_replace("Isoseq","Iso-seq"),
         CPM = CPM %>% round(digits = 2)) %>% arrange(TE)

merge_chimera_order2_order <- All_TE_merge2 %>% select(TE, type, my_col) %>% distinct()
All_TE_merge2_FC <- All_TE_merge2 %>% 
  pivot_wider(names_from = "KD", values_from = "CPM") %>% 
  mutate(CAGE_FC = log2(`CAGE Piwi KD` / `CAGE EGFP KD`),
         Isoseq_FC = log2(`Iso-seq Piwi KD` / `Iso-seq EGFP KD`),
         KD = "log2FC")

draw_merge_chimera_table <- function(){
  g1 <- All_TE_merge2 %>% filter(str_detect(KD,"CAGE")) %>% mutate(KD = KD %>% str_remove("CAGE ")) %>% 
    ggplot(aes(x = KD, y = TE %>% fct_rev, fill =CPM)) +
    geom_tile()+ geom_text(aes(label = CPM))+
    scale_fill_gradient(low = "white", high = "slateblue3",trans = "log10", na.value = "gray74",
                        breaks = c(1,10,100,1000,7000), labels = c(1,10,100,1000,7000))+
    theme_minimal()+
    theme(axis.title.y = element_blank(), axis.title.x = element_blank(),legend.position = "bottom",
          legend.key.width = unit(0.6,"cm"),legend.title = element_blank(),
          axis.text.y = element_text(color = merge_chimera_order2_order$my_col %>% rev(), size = 8))

  g2 <- All_TE_merge2_FC %>% 
    ggplot() + geom_tile(aes(x = KD, y = TE %>% fct_rev, fill =CAGE_FC))+
    geom_text(aes(x = KD, y = TE %>% fct_rev, label = if_else(is.na(CAGE_FC),"NA",if_else(is.infinite(CAGE_FC),"NA",""))))+
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral") %>% rev(),
                         na.value = "gray74",breaks=c(-5,0,5),limits = c(-5,5),oob=scales::squish)+
    theme_minimal()+ 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(),
          strip.text.x = element_blank(),legend.title = element_blank(),
          legend.position = "bottom", strip.background=element_blank(), legend.key.width = unit(0.3,"cm"))
  
  g3 <- All_TE_merge2 %>% filter(str_detect(KD,"Iso-seq")) %>% mutate(KD = KD %>% str_remove("Iso-seq ")) %>% 
    ggplot(aes(x = KD, y = TE %>% fct_rev, fill =CPM)) +
    geom_tile()+ geom_text(aes(label = CPM))+
    scale_fill_gradient(low = "white", high = "slateblue3",trans = "log10", na.value = "gray74",
                        breaks = c(1,10,100,1000,7000), labels = c(1,10,100,1000,7000))+
    theme_minimal()+
    theme(axis.title.y = element_blank(), axis.title.x = element_blank(),legend.position = "bottom",
          legend.key.width = unit(0.6,"cm"),legend.title = element_blank(),
          axis.text.y = element_blank())
  
  g4 <- All_TE_merge2_FC %>% 
    ggplot() + geom_tile(aes(x = KD, y = TE %>% fct_rev, fill =Isoseq_FC))+
    geom_text(aes(x = KD, y = TE %>% fct_rev, label = if_else(is.na(Isoseq_FC),"NA",if_else(is.infinite(Isoseq_FC),"NA",""))))+
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral") %>% rev(),
                         na.value = "gray74",breaks=c(-5,0,5),limits = c(-5,5),oob=scales::squish)+
    theme_minimal()+ 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(),
          strip.text.x = element_blank(),legend.title = element_blank(),
          legend.position = "bottom", strip.background=element_blank(), legend.key.width = unit(0.3,"cm"))

  graph <- cowplot::plot_grid(g1, g2,g3,g4, align = "h", nrow = 1,axis = "tb",
                              rel_widths = c(11/24,1/8,7/24,1/8))
  png(paste0(outdir,"merged_fusion_tile.png"), width = 6.5*ppi, height = 5.5*ppi, res = ppi)
  print(graph)
  dev.off()
}

draw_merge_chimera_table()

#######################################################################################################
# CLIP-seq count
#######################################################################################################
read_info_CLIP <- read_tsv("/media/pericles/CLASH/CLIPseq/read_info.tsv") %>% 
  mutate(type = file %>% str_remove("_f.fastq.gz") %>% str_remove("Discas/"),
         num_seqs = num_seqs %>% str_remove_all(",") %>% as.double()) %>%
  select(type, num_seqs)
filedir_CLIP <- "/media/pericles/CLASH/CLIPseq/gene/"
files_CLIP <- dir(filedir_CLIP, pattern = "\\_count.txt$", full.names = TRUE) %>%
  .[1:3]
files2_CLIP <- files_CLIP %>%
  str_remove(paste0(filedir_CLIP, "/")) %>% str_remove("\\_count.txt")
#Using Regular expression, you can get short name of the file.
df <- vector("list", length(files_CLIP))
for (i in seq_along(files_CLIP)) {
  df[[i]] <- read_tsv(files_CLIP[[i]],
                      col_names = c("gene", "length", "count", "dis")) %>%
    filter(length !=0) %>% select(-contains("dis"))
  df[[i]][["type"]] <- files2_CLIP[[i]]
}
tx2_CLIP <- bind_rows(df) %>% 
  left_join(read_info_CLIP, by = "type") %>% 
  mutate(CPM = count*1000000/num_seqs) %>% select(-count, -num_seqs) %>%
  spread(key = type, value = CPM)
tx3_CLIP <- tx2_CLIP %>% separate(gene, into = c("target", "targetname", "targettype"), sep = "_")

tx3_CLIP %>% write_csv(paste0(outdir, "CLIP_Matrix.csv"))

#######################################################################################################
# genic piRNA classification, create piRNA CPM table!
#######################################################################################################
pioutdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/piRNA/"
tx2 <- read_tsv(paste0("/media/pericles/CLASH/piRNA/","piRNA_CPM.tsv"))

piRNA_tsv <- tx2 %>% select(name,seq,pilen)


files4 <- c("piRNA1","piRNA2","piRNA3","SRR2749801","SRR2749802","SRR9158321")
df2 <- vector("list", length(files4))
for (i in seq_along(files4)) {
  df2[[i]] <- piRNA_tsv %>% 
    left_join(read_tsv(paste0("/media/pericles/CLASH/piRNA/", "otherdata/piRNA_", files4[[i]], "_for_graph.tsv")) %>% 
                select(seq, CPM), by = "seq")
  df2[[i]][["sample"]] <- files4[[i]]
}

merge_CPM <- df2 %>% bind_rows() %>%
  mutate(CPM = if_else(is.na(CPM), 0, CPM)) %>%
  select(name, seq, CPM, sample, pilen) %>% 
  pivot_wider(names_from = sample, values_from = CPM,values_fill =list(CPM=0)) %>% 
  mutate(CPM_avg = (SRR2749801 + SRR2749802 + SRR9158321 + piRNA1 + piRNA2 + piRNA3) / 6)


piRNAfromorder <- c("rmst","flam","20A","TE","genic","others")
piRNAgroup <- read_tsv(paste0("/media/pericles/CLASH/piRNA/", "anno/piRNA_all_annotation.tsv"))
piRNAremove <- read_tsv("/media/pericles/CLASH/piRNA/genome/piRNA_all_remove.bed",
                        col_names = c("chr1","start1","end1","name1","dis1","dis2",
                                      "chr2","start2","end2","name2","dis3","dis4")) %>%
  select(name1) %>% distinct() %>% 
  mutate(piRNA = name1 %>% str_extract("^piRNA\\d+"),
         type = "TEmapped") %>% select(-name1)

tx4 <- merge_CPM %>% left_join(piRNAgroup, by = "name") %>% 
  left_join(piRNAremove %>% mutate(type2 = type,name = piRNA) %>% select(name,type2), by = "name") %>% 
  mutate(type3 = case_when((type =="piRNA" & is.na(type2)) ~ "others",
                           (!(type %in% c("rmstpiRNA","flampiRNA","20ApiRNA")) & type2 == "TEmapped") ~ "TE",
                           TRUE ~ type)) %>% 
  mutate(type3 = type3 %>% str_remove("piRNA") %>% str_replace("TE2","TE")) 

tx4 %>% write_tsv(paste0(pioutdir, "piRNA_CPM.tsv"))
tx4 %>% select(name, type3) %>% 
  write_tsv(paste0(pioutdir, "piRNA_classification.tsv"))


#######################################################################################################
# Isoseq count generation
#######################################################################################################
cluster_report <- read_csv("/media/ariel/Isoseq/process/polish.cluster_report.csv")

 cluster_report2 <- cluster_report %>% 
  mutate(sample = case_when(str_detect(read_id,"m54229_201205_005007") ~ "siEGFP",
                            str_detect(read_id,"m54229_201205_211937") ~ "siPIWI",
                            str_detect(read_id,"m54209_221109_230936") ~ "siEGFP2",
                            str_detect(read_id,"m54209_221110_193755") ~ "siPIWI2",
                            TRUE ~ "none")) %>% 
  group_by(cluster_id, sample) %>% summarise(count = n()) %>% ungroup()

cluster_report2 %>% write_tsv("/media/ariel/Isoseq/process/polish_cluster_count.tsv")


#######################################################################################################
# CAGE library size estimation
#######################################################################################################
ins_len_file <- dir("/media/hermione/CAGE/stat/", pattern = "_len.tsv$", full.names = TRUE)
ins_len_file2 <- ins_len_file  %>%
  str_remove("_len.tsv$") %>%
  str_remove(paste0("/media/hermione/CAGE/stat/", "/"))
lendf <- vector("list", length(ins_len_file))
for (i in seq_along(ins_len_file)) {
  lendf[[i]] <- read_tsv(ins_len_file[[i]],
                         col_names = c("len","all_pair","inner_pair","outer_pair","other")) %>% 
    filter(len < 5000) %>% filter(len != 0)
  lendf[[i]][["sample"]] <- ins_len_file2[[i]]
}
ins_len <- lendf %>% bind_rows() %>% 
  select(len, inner_pair,sample) %>% 
  mutate(tes_rep = sample %>% str_extract("_\\d") %>% str_remove("_") %>% as.integer(),
         sample = sample %>% str_remove("_\\d")) %>% 
  group_by(len, inner_pair,sample) %>%
  summarise(count = sum(inner_pair)) %>% ungroup()


ins_len %>% filter(len > 100) %>% 
  mutate(len = 5*(len%/%5) ) %>% group_by(len) %>% summarise(count = sum(count)) %>% ungroup() %>% 
  ggplot(aes(x = len, y = count))+
  geom_point()+ geom_smooth(method = "loess", span = 0.03)+
  scale_x_log10(breaks =c(100,200,300,400,500,750,1000,1500,2000,5000) )+
  labs(x = "estimated read length (nt)", y = "count")+ theme_minimal() + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.grid.minor = element_blank())+
  ggsave(paste0(outdir,"estimated_read_length_all.png"), width =3.5, height =2.8, dpi = 600)+
  ggsave(paste0(outdir,"estimated_read_length_all.pdf"), width =3.5, height =2.8)

#######################################################################################################
# CAGE QC2
#######################################################################################################
# TSS startsite

library(ggseqlogo)
motifdir <- "/media/hermione/CAGE/QC/motif/"
ctssdir <- "/media/hermione/CAGE/ctss/"
file_cage <- dir(motifdir, pattern = ".tsv$", full.names = TRUE)
file_cage_ctss <- dir(motifdir, pattern = ".bed$", full.names = TRUE)
file_cage2 <- file_cage  %>%
  str_replace_all(".tsv$", "") %>%
  str_replace(paste0(motifdir, "/"), "")
cagedf <- vector("list", length(file_cage))
for (i in seq_along(file_cage)) {
  cagedf[[i]] <- read_tsv(file_cage[[i]],
                          col_names = c("name","seq"))
  cagedf[[i]][["sample"]] <- file_cage2[[i]]
}
cagedf2 <- vector("list", length(file_cage_ctss))
for (i in seq_along(file_cage_ctss)) {
  cagedf2[[i]] <- read_tsv(file_cage_ctss[[i]],
                           col_names = c("chr","start","end","name","cov","strand"))
  cagedf2[[i]][["sample"]] <- file_cage2[[i]]
}
cage_seq <- cagedf %>% bind_rows() %>%
  separate(name, into = c("chr","name"),sep = ":") %>%
  separate(name, into = c("name","strand"),sep = "\\(") %>%
  separate(name, into = c("start","end"),sep = "-") %>%
  mutate(start = start %>% as.integer(),end = end %>% as.integer(), strand = strand %>% str_remove("\\)")) %>%
  mutate(seq = if_else(strand == "-" ,str_sub(seq,start = 1, end = -2),str_sub(seq,start = 2, end = -1)))
cage_all <- cagedf2 %>% bind_rows()

cage_all2 <- cage_all %>%
  left_join(cage_seq, by = c("chr","start","end","strand","sample"))

cage_all3 <- cage_all2 %>% group_by(sample) %>%
  mutate(cov2 = cov * 1000000 /sum(cov),
         intensity = case_when(cov2>5 ~ "strong", cov2 > 1 ~ "medium",TRUE ~ "weak")) %>%
  ungroup()

logo <- cage_all3 %>% filter(cov >=5) %>% filter(sample == "siEGFP_rep1") %>%
  mutate(intensity = paste0(intensity," TSS")) %>%
  select(seq,intensity) %>%
  distinct() %>%
  mutate(seq2 = seq %>% str_sub(start = 1L, end = 8L),
         type2 = intensity) %>%
  select(type2, seq2)
logo2 <- logo %>% group_by(type2) %>% mutate(total = n()) %>% ungroup() %>%
  mutate(type2 = paste0(type2, " (n=",total, ")")) %>% select(-total)

tx <- logo2 %>% group_nest(type2) %>% 
  mutate(order = case_when(str_detect(type2, "strong") ~ 1,
                           str_detect(type2, "medium") ~ 2,
                           str_detect(type2, "weak") ~ 3,
                           TRUE ~4)) %>% 
  arrange(order) %>% select(-order)
tx2 <- map(1:length(tx$type2), function(x) pull(tx$data[[x]]))
names(tx2) <- tx$type2
ggseqlogo(tx2,ncol = 1)+
  scale_x_continuous(breaks=1:8, label = c("-4","-3","-2","-1","1","2","3","4"))+
  ggsave(paste0(outdir,"base_motif_logo_all.png"), width =2.8, height =4.8, dpi = 600)+
  ggsave(paste0(outdir,"base_motif_logo_all.pdf"), width =2.8, height =4.8)

#######################################################################################################
# piRNA abundance
#######################################################################################################
library(stringi)
library(cowplot)

pifigoutdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/piRNA/abundance/"
piRNA_CPM_cal_d <- read_tsv("/media/pericles/CLASH/piRNA/read_info.tsv") %>% 
  mutate(name = file %>% str_remove(".fa")) %>% 
  select(name, num_seqs) %>%  filter(name == "piRNA3")

TE_anno <- read_tsv("/media/pericles/TEfind/database/annotateTE/TE_bed6.bed",
                    col_names = c("target","start","end","des","dis","direction")) %>%
  mutate(target = target %>% str_replace("-", "\\|") %>% str_remove("\\-"),
         des = des %>% str_remove("CDS_") %>%
           str_replace("five_prime_","5'") %>% str_replace("three_prime_","3'"))

TE <- read_tsv("/media/pericles/TEfind/database/TE_dm.tsv", col_names = c("TE3","fasta")) %>%
  mutate(target = TE3 %>% str_replace("\\-","|") %>% str_remove("\\-"),
         totallen=str_length(fasta),
         start_TE=0) %>%
  select(target,start_TE,totallen)

piRNA_variety <- read_tsv("/media/pericles/CLASH/piRNA/TE/mismatch/piRNA_all_3_reverse_1.bedgraph",
                          col_names = c("chr","start","end","cov"))

piRNA_express <- read_tsv("/media/pericles/CLASH/piRNA/TE/mismatch/piRNA3_3_reverse_1.bedgraph",
                          col_names = c("chr","start","end","cov")) %>% 
  mutate(name = "piRNA3") %>% left_join(piRNA_CPM_cal_d, by = "name") %>% 
  mutate(CPM = 1000000*cov/num_seqs)

flam_TE1 <- read_tsv("/media/pericles/TEfind/flam/rmblast/flam_hap1_TE_minus.bed",
                     col_names = c("chr", "start","end","name","score","strand")) %>% 
  mutate(sample = "flam hap1")
flam_TE2 <- read_tsv("/media/pericles/TEfind/flam/rmblast/flam_hap2_TE_minus.bed",
                     col_names = c("chr", "start","end","name","score","strand")) %>% 
  mutate(sample = "flam hap2")
l20A_TE1 <- read_tsv("/media/pericles/TEfind/flam/rmblast/20A_hap1_TE_minus.bed",
                     col_names = c("chr", "start","end","name","score","strand")) %>% 
  mutate(sample = "20A hap1")
l20A_TE2 <- read_tsv("/media/pericles/TEfind/flam/rmblast/20A_hap2_TE_minus.bed",
                     col_names = c("chr", "start","end","name","score","strand"))%>% 
  mutate(sample = "20A hap2")
flam_TE2 <- bind_rows(flam_TE1,flam_TE2,l20A_TE1,l20A_TE2 ) %>% 
  separate(name, sep = "_" ,into = c("name", "TEstart","TEend")) %>% 
  mutate(target = name %>% str_replace("\\-","|") %>% str_remove("\\-"),
         TEstart = TEstart %>% as.integer(),
         TEend = TEend %>% as.integer())

piRNA_expression_flam <- function(data, fil,out,height_in) {
  g4 <- ggplot() +
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=-0.5,ymax=5.5),fill = "transparent")+
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=0,ymax=5),fill = "darkgray")+
    geom_rect(data = TE_anno %>% filter(target==fil),aes(xmin=start,xmax = end,ymin=0,ymax=5,fill = des),color = "bisque4")+
  #  geom_text(data = TE_anno %>% filter(target==fil),aes(x=(start+end)/2,y=2.5,label=des))+
    geom_text_repel(data = subset(TE_anno, target==fil),
                    aes(label = des,x=(start+end)/2,y=2.5), color = "black", size = 3, min.segment.length = 0,
                    segment.color = "gray90", nudge_x = 0,nudge_y = 0,point.size = NA)+
    theme_minimal()+ labs(x = out, y="",fill="",color = "")+scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))+
    guides(fill="none") +theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(), panel.grid.major.x  = element_line(color = "black"))
  g3 <- flam_TE2 %>% filter(target==fil) %>% mutate(num = row_number() %>% as.integer()) %>% 
    ggplot() +
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=0,ymax=1),fill = "transparent")+
    geom_rect(aes(xmin=TEstart,xmax = TEend,ymin=num-0.95,ymax=num-0.05, fill = sample, color = sample)) +
    theme_minimal()+ labs(y="flam/20A insertion \n(antisense)",fill="",color = "") +scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))+
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(),panel.grid.minor.y = element_blank())
  
  g1 <- piRNA_express %>% 
    mutate(target = chr %>% str_replace("\\-","|") %>% str_remove("\\-")) %>% 
    filter(target==fil) %>%
    ggplot()+
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=0,ymax=1),fill = "transparent")+
    geom_step(aes(x=start,y=-CPM),color = "blue")+
    labs(y="piRNA expression (CPM)",fill="",color = "")+
    theme_minimal()+scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))+
    theme(legend.position = c(0.8,0.8),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  g2 <- piRNA_variety %>% 
    mutate(target = chr %>% str_replace("\\-","|") %>% str_remove("\\-")) %>% 
    filter(target==fil) %>%
    ggplot()+
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=0,ymax=1),fill = "transparent")+
    geom_step(aes(x=start,y=-cov),color = "red")+
    labs(y="piRNA variety (count)",fill="",color = "")+
    theme_minimal()+scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))+
    theme(legend.position = c(0.8,0.8),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  graph <- plot_grid(g1,g2,g3,g4,  align = "v",axis = "lr", nrow = 4, rel_heights = height_in)
  png(paste0(pifigoutdir,"piRNAcov_",out,".png"), width = 4.5*ppi, height = 4.5*ppi, res = ppi)
  print(graph)
  dev.off()
}

piRNA_expression_flam(fil="Dmel|blood",out = "blood",height_in=c(1/3,1/3,2/15, 1/5))
piRNA_expression_flam(fil="Dmel|gypsy",out = "gypsy",height_in=c(1/4,1/4,3/10,1/5))
piRNA_expression_flam(fil="Dmel|mdg1",out = "mdg1",height_in=c(1/4,1/4,3/10,1/5))
piRNA_expression_flam(fil="Dmel|ZAM",out = "ZAM",height_in=c(1/3,1/3,2/15, 1/5))


# mdg1_zoom

piRNA_expression_flam2 <- function(data, fil,out,height_in,zoom_in1,zoom_in2) {
  g4 <- ggplot() +
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=-0.5,ymax=5.5),fill = "transparent")+
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=0,ymax=5),fill = "darkgray")+
    geom_rect(data = TE_anno %>% filter(target==fil),aes(xmin=start,xmax = end,ymin=0,ymax=5,fill = des),color = "bisque4")+
    #  geom_text(data = TE_anno %>% filter(target==fil),aes(x=(start+end)/2,y=2.5,label=des))+
    geom_text_repel(data = subset(TE_anno, target==fil),
                    aes(label = des,x=(start+end)/2,y=2.5), color = "black", size = 3, min.segment.length = 0,
                    segment.color = "gray90", nudge_x = 0,nudge_y = 0,point.size = NA)+
    theme_minimal()+ labs(x = out, y="",fill="",color = "")+scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))+
    guides(fill="none") +theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(), panel.grid.major.x  = element_line(color = "black"))
  g3 <- flam_TE2 %>% filter(target==fil) %>% mutate(num = row_number() %>% as.integer()) %>% 
    ggplot() +
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=0,ymax=1),fill = "transparent")+
    geom_rect(aes(xmin=TEstart,xmax = TEend,ymin=num-0.95,ymax=num-0.05, fill = sample), color= "bisque4") +
    theme_minimal()+ labs(y="flam/20A insertion \n(antisense)",fill="",color = "") +scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))+
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(),panel.grid.minor.y = element_blank())
  
  g1 <- piRNA_express %>% 
    mutate(target = chr %>% str_replace("\\-","|") %>% str_remove("\\-")) %>% 
    filter(target==fil) %>%
    ggplot()+
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=0,ymax=1),fill = "transparent")+
    geom_step(aes(x=start,y=-CPM),color = "blue")+
    labs(y="piRNA expression (CPM)",fill="",color = "")+  coord_cartesian(ylim = c(0,zoom_in1))+
    theme_minimal()+scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))+
    theme(legend.position = c(0.8,0.8),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  g2 <- piRNA_variety %>% 
    mutate(target = chr %>% str_replace("\\-","|") %>% str_remove("\\-")) %>% 
    filter(target==fil) %>%
    ggplot()+
    geom_rect(data = TE %>% filter(target==fil), aes(xmin=start_TE,xmax = totallen,ymin=0,ymax=1),fill = "transparent")+
    geom_step(aes(x=start,y=-cov),color = "red")+
    labs(y="piRNA variety (count)",fill="",color = "")+ coord_cartesian(ylim = c(0,zoom_in2))+
    theme_minimal()+scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))+
    theme(legend.position = c(0.8,0.8),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  graph <- plot_grid(g1,g2,g3,g4,  align = "v",axis = "lr", nrow = 4, rel_heights = height_in)
  png(paste0(pifigoutdir,"piRNAcov_",out,"_zoom.png"), width = 4.5*ppi, height = 4.5*ppi, res = ppi)
  print(graph)
  dev.off()
}


piRNA_expression_flam2(fil="Dmel|mdg1",out = "mdg1",height_in=c(1/4,1/4,3/10,1/5), zoom_in1 = 600, zoom_in2 = 150)




