
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/CLASH/vienna/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gplots)))
suppressMessages(suppressWarnings(require(RColorBrewer)))
suppressMessages(suppressWarnings(require(gtools)))
suppressMessages(suppressWarnings(require(gridExtra)))
suppressMessages(suppressWarnings(require(stringi)))
suppressMessages(suppressWarnings(require(gt)))
suppressMessages(suppressWarnings(require(ggpubr)))

ppi <- 600
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

outdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/fig2/"
CLASHdir <- "/media/pericles/CLASH/"
CLASH_shufdir <- "/media/pericles/CLASH/shuffle/"

#######################################################################################################
# clustering
# fig2A&B&D&sup_figA&E
#######################################################################################################
referenceRNAorder <- c("TE","mRNA","ncRNA","rRNA","snoRNA","snRNA","tRNA","pseudogene")

Vienna_sle <-read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/", "CLASH_sup_file1.csv"),
                      col_types = cols(dG = col_double()))


Vienna2 <- Vienna_sle %>% 
  select(group, sample, piRNA,targettype, pilen,piRNAseq, piRNAmatch, dG, readname) %>%
  filter(is.finite(dG)) %>%
  mutate(heat = piRNAmatch %>% str_replace_all("\\.","0-") %>% str_replace_all("\\(","1-")) %>%
  mutate(score = map(heat, function(x) {x %>% str_sub(end = -2) %>% str_split(pattern = "-")}),
         index2 = row_number()) %>%
  unnest(cols = score) %>%
  mutate(base = map(score, function(x) {seq(1, length(x), by = 1)})) %>%
  select(group, sample,index2, dG, pilen, piRNA, score, base,targettype, readname) %>%
  unnest(cols = c(score, base)) %>%
  mutate(score_G = as.integer(score) * dG,
         gene_name = paste0(piRNA, "_", index2))

kmeans_cal <- function(data, center){
  kmeans <- Vienna2 %>%
    select(-score, -piRNA, -group,-sample, -targettype) %>%
    filter(pilen %in% c(23,24,25, 26, 27, 28,29)) %>%
    group_by(pilen) %>% nest() %>%
    mutate(miranda = map(data, function(x) {x %>% spread(base, score_G)}),
           km = map(miranda, function(x) {x %>% select(-index2, -gene_name, -dG, -readname) %>% scale() %>%
               kmeans(centers = center, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")}),
           cluster = map(km, function(x) {x[["cluster"]]})) %>%
    select(-km, -data) %>%
    unnest(cols = c(miranda, cluster))
  kmeans2 <- kmeans %>% ungroup()%>% select(gene_name, cluster) %>%
    inner_join(Vienna2, by = "gene_name") %>%
    group_by(pilen,cluster) %>% nest() %>%
    mutate(center_dG = map(data, function(x) {x$dG %>% mean()})) %>%
    unnest(cols = c(data, center_dG)) %>%
    arrange(pilen, cluster, gene_name, base) %>%
    mutate(piRNAlen2 = paste0("piRNA (", pilen, "nt)"))
  kmeans3 <- kmeans2 %>% select(pilen, cluster, center_dG) %>%
    distinct() %>%
    group_by(pilen) %>% nest_legacy() %>%
    mutate(data2 = map(data, function(x) {x %>% arrange(center_dG) %>%
        mutate(cluster2=row_number() %>% as.integer())})) %>% select(-data) %>%
    unnest_legacy() %>%
    select(-center_dG)
  kmeans4 <- kmeans2 %>% left_join(kmeans3, by = c("pilen", "cluster")) %>%
    ungroup() %>%
    mutate(cluster = paste0("cluster", cluster2)) %>%
    select(-cluster2)
  kmeans4 %>% write_csv(paste0(outdir, "kmeans_center3_output.csv"))
  return(kmeans4)
}

draw_kmeans <- function(kmeans4, center){
  tiledata <- kmeans4 %>%
    select(cluster,pilen, gene_name, piRNA, group, sample,
           readname,piRNAlen2,dG,index2,targettype,center_dG) %>%
    distinct()
  kmeans4 %>% ggplot() +
    geom_raster(aes(y =fct_reorder(gene_name, center_dG) %>% fct_rev(), fill = score_G, x=base))+
    facet_wrap(~piRNAlen2, strip.position ="top", scales = "free_y", ncol=7)+
    geom_point(data = tiledata, aes(x = 0, y =gene_name, color = cluster), size = 0.5, shape = 15) +
    #guides(color = "none")+
    scale_x_continuous(breaks = c(1,10,20))+theme_minimal()+
    scale_fill_gradient2(low = "black",mid = "black", high = "white", midpoint = -50)+
    labs(fill = "\u0394G (kcal/mol)", color = "",x = "bases (nt)", y= paste0("ordered by k-means (center = ",center,")"))+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "top",
          strip.background=element_blank(), panel.grid = element_blank())+
    ggsave(paste0(outdir,"fig2A.png"), width =7, height =5, dpi = 600)+
    ggsave(paste0(outdir,"fig2A.pdf"), width =7, height =5)
  kmeans4 %>% select(pilen, gene_name, dG, cluster, group, sample, piRNAlen2,targettype) %>%
    distinct() %>%
    mutate(targettype = if_else(targettype %in% c("mRNA(TE)","ensemblTE","rmblastTE","pbsvTE"), "TE",targettype)) %>%
    group_by(piRNAlen2, cluster) %>%
    count(targettype) %>% mutate(per = 100 *n/sum(n)) %>% ungroup() %>%
    mutate(targettype = targettype %>% factor(levels = referenceRNAorder)) %>% 
    ggplot(aes(x = cluster %>% fct_rev(), y = per, fill = targettype)) +
    geom_col(color = "black")+
    theme_minimal() + coord_flip() +
    facet_wrap(~piRNAlen2, ncol=7) + 
    labs(y = "Fraction (%)") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          axis.title.y = element_blank()) +
    ggsave(paste0(outdir,"fig2B.png"), width =10, height =0.6*center+1.2, dpi = 300)+
    ggsave(paste0(outdir,"fig2B.pdf"), width =10, height =0.6*center+1.2)
  kmeans4 %>% select(pilen, gene_name, dG, cluster, group, sample, piRNAlen2) %>%
    distinct() %>%
    ggplot(aes(x = cluster, y = dG, fill = cluster))+
    geom_boxplot()+
    facet_wrap(~piRNAlen2, ncol=7)+
    labs(fill = "", y = "dG (kcal/mol)", x = "") +
    theme_minimal()+
    theme(legend.position = "none") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    ggsave(paste0(outdir,"sup_fig2E.png"), width =3*center, height =3, dpi = 300)+
    ggsave(paste0(outdir,"sup_fig2E.pdf"), width =3*center, height =3)
}


Vienna2 %>% kmeans_cal(center = 3) %>% draw_kmeans(center = 3)


elbow <- function(data) {
  kmeans_elbow <- data %>%
    filter(pilen == 26) %>%
    select(-score, -piRNA, -group,-sample, -targettype, -pilen, -readname) %>%
    spread(base, score_G) %>%
    select(-index2, -dG, -gene_name)
  pct_var <- data.frame(pct_var = 0, num_clusters = 2:20)
  totalss <- kmeans(kmeans_elbow, centers = 20, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")$totss
  for(i in 2:20){
    pct_var[i-1, 'pct_var'] <- kmeans(kmeans_elbow, centers = i, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")$betweenss/totalss
  }
  return(pct_var)
}
elbow_result <- Vienna2 %>% elbow()

ggplot(elbow_result, aes(x = num_clusters, y = pct_var*100)) +
  geom_line() + geom_point() + theme_minimal() +
  ylim(0, 100) + xlab("Number of clusters") + ylab("% Variance Explained") +
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  scale_x_continuous(breaks = c(2,3,4,5,10,15,20)) +
  ggsave(paste0(outdir,"sup_fig2A.png"),width =2.5, height =2, dpi = 600)+
  ggsave(paste0(outdir,"sup_fig2A.pdf"),width =2.5, height =2)

elbow_result2 <- elbow_result %>% 
  mutate(pct_var2 = append(0,head(pct_var, -1)), dif_pct = pct_var - pct_var2) %>% 
  filter(num_clusters < 11) %>%filter(num_clusters >2) %>% 
  mutate(samplename = paste0((num_clusters -1), " -> ",num_clusters))

elbow_result2 %>% 
ggplot(aes(x = 100 * dif_pct, y = samplename)) +
  geom_col(fill = "deepskyblue4",color = "black") + theme_minimal() +
  xlab("Gain in Variance explained (%)") + ylab("Number of clusters") +
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
  ggsave(paste0(outdir,"sup_fig2A_2.png"),width =3.5, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"sup_fig2A_2.pdf"),width =3.5, height =2.5)


# read k=3 file
kmeans4 <- read_csv(paste0(outdir, "kmeans_center3_output.csv"), col_types = cols(dG = col_double()))

kmeans5 <- kmeans4 %>% select(piRNAlen2, cluster,group,sample,gene_name, dG, score) %>%
  group_by(group, sample, cluster,dG,piRNAlen2, gene_name) %>%
  nest_legacy() %>%
  mutate(score = map(data, function(x) {x$score %>% paste0(collapse = "_")})) %>%
  select(-data) %>% unnest_legacy() %>%
  mutate(mismatch = str_count(score, "0"))

kmeans5 %>% group_by(cluster, mismatch) %>% summarise(count = n()) %>% ungroup() %>%
  ggplot(aes(x = mismatch, y = count, color = cluster))+
  geom_line()+geom_area(aes(fill = cluster),alpha = 0.1, position = "identity", color = NA) +
  xlab("mismatch count") + ylab("Number of chimeric reads") +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill = "none") +
  theme_minimal() +
  theme(legend.position = c(.8,.8), legend.title = element_blank())+
  ggsave(paste0(outdir,"fig2D.png"), width =3, height =3, dpi = 300)+
  ggsave(paste0(outdir,"fig2D.pdf"), width =3, height =3)

# create all containing files....

tiledata <- kmeans4 %>%
  select(cluster,pilen, gene_name, piRNA, group, sample,
         readname,piRNAlen2,dG,index2,targettype,center_dG) %>%
  distinct() %>% arrange(piRNAlen2)
kmeans4 %>% ggplot() +
  geom_raster(aes(y =fct_reorder(gene_name, pilen) , fill = score_G, x=base))+
  geom_point(data = tiledata, aes(x = 0, y =gene_name, color = cluster), size = 0.5, shape = 15) +
  facet_wrap(~cluster, strip.position ="top",scales = "free_y", ncol=3)+
  guides(color = "none")+
  scale_x_continuous(breaks = c(1,10,20))+ theme_minimal()+
  scale_fill_gradient2(low = "black",mid = "black", high = "white", midpoint = -50)+
  labs(fill = "\u0394G (kcal/mol)", color = "",x = "bases (nt)", y= paste0("ordered by k-means (center = ",3,")"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        strip.background=element_blank(), panel.grid = element_blank())+
  ggsave(paste0(outdir,"fig2A_all.png"), width =5, height =4, dpi = 600)

kmeans4 %>% select(pilen, gene_name, dG, cluster, group, sample, piRNAlen2,targettype,readname) %>%
  distinct() %>%
  mutate(readcount = readname %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer()) %>% 
  mutate(targettype = if_else(targettype %in% c("mRNA(TE)","ensemblTE","rmblastTE","pbsvTE"), "TE",targettype)) %>%
  group_by(cluster,targettype) %>%
  summarise(n = sum(readcount)) %>% ungroup() %>% group_by(cluster) %>% mutate(prop = 100 *n/sum(n)) %>% 
  mutate(ypos = cumsum(prop)-0.5*prop) %>% ungroup() %>% 
  mutate(targettype = targettype %>% factor(levels = referenceRNAorder %>% fct_rev())) %>% 
  ggplot(aes(x = "", y = prop, fill = targettype)) +
  geom_bar(stat="identity", width=1, color="white", linewidth=0.1) +
  coord_polar("y", start=0, direction = 1) +
  theme_void() +
  facet_wrap(~cluster, ncol=3) + 
  guides(fill= guide_legend(title = NULL))+
  geom_text(aes(label = if_else(prop > 3,paste0(round(prop,1), "%"),"")), position = position_stack(vjust = 0.5), size= 3)+
  theme(legend.position = "bottom") +
  ggsave(paste0(outdir,"fig2B_pie_all.png"), width =5, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"fig2B_pie_all.pdf"), width =5, height =2.5)

RNAtable_CLASH_cluster <- kmeans4 %>% select(pilen, gene_name, dG, cluster, group, sample, piRNAlen2,targettype,readname) %>%
  distinct() %>%
  mutate(readcount = readname %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer()) %>% 
  mutate(targettype = if_else(targettype %in% c("mRNA(TE)","ensemblTE","rmblastTE","pbsvTE"), "TE",targettype)) %>%
  group_by(cluster,targettype) %>%
  summarise(n = sum(readcount)) %>% ungroup() %>% group_by(cluster) %>% mutate(prop = 100 *n/sum(n)) %>% 
  mutate(ypos = cumsum(prop)-0.5*prop) %>% ungroup()
RNAtable_CLASH_cluster %>% 
  write_csv(paste0(outdir,"CLASH_RNAtype_cluster_table.csv"))

#######################################################################################################
# fig1C
# comparing shuffle
#######################################################################################################
Vienna_shuffle <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/", "CLASH_sup_file3.csv"),
                           col_types = cols(dG = col_double())) %>% 
  mutate(dG = if_else(!(str_detect(piRNAmatch, "\\(")), Inf, dG)) %>%
  select(sample, readname,dG,piRNAmatch,targetmatch) %>%
  rename(dG_shuf = dG ,
         piRNAmatch_shuf = piRNAmatch,
         targetmatch_shuf = targetmatch) %>%
  right_join(Vienna_sle, by = c("sample","readname"))

Vienna_shuffle2 <- Vienna_shuffle %>% select(sample,readname,dG,piRNAmatch) %>%
  mutate(type = "chimera") %>%
  bind_rows(Vienna_shuffle %>% select(sample,readname,dG_shuf,piRNAmatch_shuf) %>%
              rename(dG =dG_shuf,piRNAmatch=piRNAmatch_shuf) %>%
              mutate(type = "shuffle")) %>%
  rename(match = piRNAmatch)

Vienna_shuffle4 <- Vienna_shuffle2 %>% mutate(pilen=str_length(match)) %>%
  mutate(heat = match %>% str_replace_all("\\.","0-") %>% str_replace_all("\\(","1-") %>% str_replace_all("\\)","1-"),
         dG = if_else(is.infinite(dG), 0,dG)) %>%
  mutate(score = map(heat, function(x) {x %>% str_sub(end = -2) %>% str_split(pattern = "-")}),
         index2 = row_number()) %>%
  unnest(cols = score) %>%
  mutate(base = map(score, function(x) {seq(1, length(x), by = 1)})) %>%
  select(sample,index2, dG,score, base,readname,type,  pilen) %>%
  unnest_legacy()

Vienna_shuffle4 %>% left_join(kmeans4 %>% select(readname, pilen, sample,cluster) %>% distinct(), by =c("readname","sample","pilen")) %>%
  filter(!(is.na(cluster))) %>%
  filter(pilen %in% c(23,24,25,26,27,28,29)) %>%
  mutate(score = score %>% as.double(),
         pilen = paste0(pilen, "nt"),
         type = if_else(type == "chimera", "CLASH-chimeras", "shuffled-chimeras")) %>%
  group_by(type,pilen, base,cluster) %>%
  summarise(base_pair = mean(score)) %>% ungroup() %>%
  ggplot(aes(x = base, y = 100 *base_pair, color =pilen)) +
  geom_line(aes(linetype = type)) + annotate("rect", xmin = 2,xmax = 8,ymin=0,ymax = 100,alpha = 0.1, fill = "blue")+
  theme_minimal()+ guides(linetype=guide_legend(nrow = 2,byrow = TRUE, keyheight = 0.9)) +
  facet_wrap(~cluster, nrow = 1)+
  coord_cartesian(xlim = c(1, 29), ylim = c(0,100)) +
  scale_x_continuous(breaks = c(1,5,10,15,20,25,29)) +
  labs(color = "", x = "bases (nt)", y= "base pairing (%)", linetype ="") +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank())+
  ggsave(paste0(outdir,"fig2C.png"),  width =5, height =3.2, dpi = 300)+
  ggsave(paste0(outdir,"fig2C.pdf"),  width =5, height =3.2)

#######################################################################################################
# fig2F
#######################################################################################################
Vienna_sle <-read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/", "CLASH_sup_file1.csv"),
                      col_types = cols(dG = col_double()))
kmeans4 <- read_csv(paste0(outdir, "kmeans_center3_output.csv"), col_types = cols(dG = col_double()))
cage_matrix <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "CAGE_Matrix.csv"))

Vienna_select <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/","CLASH_sup_file2.csv"),
                          col_types = cols(dG = col_double()))

kmeans_cluster <- kmeans4 %>% select(sample,readname, cluster) %>% distinct()
Vienna_cluster <- Vienna_select %>% left_join(kmeans_cluster, by = c("sample","readname")) %>%
  select(sample, group, cluster, readname, target, targetname, targettype, bedstart, bedend, piRNA, piRNAtype,dG,piRNAmatch,targetmatch)
Vienna_cluster %>% write_csv(paste0(outdir, "clash_selection_cluster.csv"))
# start here
Vienna_cluster <- read_csv(paste0(outdir, "clash_selection_cluster.csv"), col_types = cols(dG = col_double()))
Vienna_cluster2 <- Vienna_cluster %>% filter(targettype!="mRNA(TE)") %>%
  select(sample,readname, cluster,target,dG) %>%
  mutate(cluster = if_else(is.na(cluster), "non",cluster)) %>%
  group_by(target) %>% mutate(totalcount = n()) %>% ungroup() %>% filter(totalcount >= 2) %>% 
  group_by(target, cluster) %>%
  mutate(count = n(), meandG = mean(dG)) %>% ungroup() %>% 
  group_by(target) %>% slice(which.max(count)) %>% ungroup() %>%
  select(target,cluster,dG) %>% distinct()

Vienna_cluster3 <- Vienna_cluster2 %>% 
  separate(target, sep = "%", into = c("target","dis1","dis2")) %>% filter(!(dis1 %in% c("pbsv","rmblast"))) %>% select(-contains("dis")) %>% 
  mutate(target = target %>% str_replace("element","-element") %>% str_replace("Beagle","-Beagle") %>% str_replace("TART","TART-"))

Vienna_cluster3 %>% 
  write_csv(paste0(outdir,"RNA_cluster_class.csv"))
  

# CAGE-seq
cage_edgeR <- cage_matrix %>% select(genename, sample,count) %>% 
  group_by(genename, sample) %>% summarise(count = sum(count)) %>% ungroup() %>% arrange(genename,sample) %>% 
  pivot_wider(names_from = "sample",values_from = "count" , values_fill =list(count=0)) 

CLASH_RNA <- cage_edgeR %>% mutate(target = genename %>% str_replace("-", "|") %>% str_replace("_", "|")) %>% 
  left_join(Vienna_cluster3,by="target") %>%
  mutate(cluster = if_else(is.na(cluster), "else",cluster))

CLASH_RNA2 <- CLASH_RNA %>% mutate(siEGFP = siEGFP_rep1 + siEGFP_rep2+siEGFP_rep3,
                                   siPIWI = siPIWI_rep1 + siPIWI_rep2+siPIWI_rep3) %>% 
  select(-target) %>% filter(cluster %in% c("cluster1","cluster2","cluster3")) %>% 
  mutate(total = rowSums(across(contains("rep"))),
         total2 = sum(total),
         total_CPM = 1000000*total/total2) %>% 
  filter(total_CPM >= 1) %>% 
  select(-total, -total2,-total_CPM) %>% mutate(delta = log2((siPIWI+1)/(siEGFP+1)))

my_comparisons <- list( c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3") )
count_all <- CLASH_RNA2  %>% 
  group_by(cluster) %>% count()

CLASH_RNA2  %>%
  ggplot(aes(x=cluster,y=delta,fill = cluster))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.8, 6.5, 5.8), method = "wilcox.test", size = 3.5 ,label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 2.2e-16, 0.0001, 0.001, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))
                    )+
 # geom_text( y = -5, data = count_all, aes(x = cluster,label = paste0("italic(n) ==",n)), parse = TRUE,size = 3) +
 # geom_text( y = -4, data = count_all, aes(x = cluster,label = paste0("n=",n)),size = 3) +
  labs(y = "log2((Piwi KD+1)/(EGFP KD+1))") +
  theme_minimal()+
  guides(fill = "none")+
  theme(axis.title.x = element_blank())+
  ggsave(paste0(outdir,"fig2F.png"),  width =2.5, height =4, dpi = 300)+
  ggsave(paste0(outdir,"fig2F.pdf"),  width =2.5, height =4)
  # n=45, n=309, n=1209
# Iso-seq
Isoseq_matrix <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "Isoseq_Matrix.csv"))
Isoseq_edgeR <- Isoseq_matrix %>% rename(genename = associated_gene) %>% select(genename, fl_EGFP, fl_PIWI)

Isoseq_RNA <- Isoseq_edgeR %>% mutate(target = genename %>% str_replace("-", "|") %>% str_replace("_", "|")) %>% 
  left_join(Vienna_cluster3,by="target") %>%
  mutate(cluster = if_else(is.na(cluster), "else",cluster))

Isoseq_RNA2 <- Isoseq_RNA %>% 
  select(-target) %>% filter(cluster %in% c("cluster1","cluster2","cluster3")) %>% 
  mutate(total = fl_EGFP + fl_PIWI,
         total2 = sum(total),
         total_CPM = 1000000*total/total2) %>% 
  filter(total_CPM >= 1) %>% 
  select(-total, -total2,-total_CPM) %>% mutate(delta = log2((fl_PIWI+1)/(fl_EGFP+1)))

my_comparisons <- list( c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3") )
Isoseq_count_all <- Isoseq_RNA2  %>% 
  group_by(cluster) %>% count()

Isoseq_RNA2  %>%
  ggplot(aes(x=cluster,y=delta,fill = cluster))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, label.y = c(8, 8.7, 8), method = "wilcox.test", size = 3.5 ,label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 2.2e-16, 0.0001, 0.001, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))
  )+
#  geom_text( y = -4.7, data = Isoseq_count_all, aes(x = cluster,label = paste0("italic(n) ==",n)), parse = TRUE,size = 3) +
#  geom_text( y = -4.7, data = Isoseq_count_all, aes(x = cluster,label = paste0("n=",n)), size = 3) +
  labs(y = "log2((Piwi KD+1)/(EGFP KD+1))") +
  theme_minimal()+
  guides(fill = "none")+
  theme(axis.title.x = element_blank())+
  ggsave(paste0(outdir,"fig2F_Isoseq.png"),  width =2.5, height =4, dpi = 300)+
  ggsave(paste0(outdir,"fig2F_Isoseq.pdf"),  width =2.5, height =4)


#######################################################################################################
# sup_fig2B
# stat of clustering
#######################################################################################################

kmeans_cluster <- kmeans4 %>% select(sample,readname, cluster) %>% distinct()
Vienna_cluster_sle <- Vienna_sle %>% left_join(kmeans_cluster, by = c("sample","readname"))
Vienna_cluster_sle %>% write_csv(paste0(outdir, "clash_sle_cluster3.csv"))



make_table <- function(data, refdata = read_count3,groupname = "sample", out) {
  data2 <- data
  groupname2 <- rlang::enquo(groupname)
  groupname3 <- rlang::enquos(groupname)
  Vienna_stat <- data2 %>% group_by(!!!groupname3) %>%
    summarise(count = n(),
              mean_pilen = mean(pilen) %>% round(2),
              sd_pilen = sd(pilen) %>% round(2),
              mean_tag = mean(targetlen) %>% round(2),
              sd_tag = sd(targetlen) %>% round(2),
              mean_mismatch = mean(str_count(piRNAmatch,"\\.")) %>% round(2),
              sd_mismatch = sd(str_count(piRNAmatch,"\\.")) %>% round(2)) %>% ungroup()
  Vienna_stat2 <- data2 %>%
    filter(is.finite(dG)) %>%
    group_by(!!!groupname3) %>%
    summarise(mean_dG = mean(dG) %>% round(2),
              sd_dG = sd(dG) %>% round(2)) %>% ungroup() %>% select(-!!groupname2)
  Vienna_stat3 <- data2 %>%
    filter(targettype=="ensemblTE") %>%
    group_by(!!!groupname3) %>%
    summarise(TEcount = n()) %>% ungroup() %>% select(-!!groupname2)
  if (dim(refdata)[1]== 0){
    Vienna_stat_all <- bind_cols(Vienna_stat,Vienna_stat2,Vienna_stat3) %>%
      mutate(anno = !!groupname2,
             pilen = paste0(mean_pilen,"±",sd_pilen),
             mismatch = paste0(mean_mismatch,"±",sd_mismatch),
             tag = paste0(mean_tag,"±",sd_tag),
             dG = paste0(mean_dG,"±",sd_dG),) %>%
      select(-contains("mean"), -contains("sd"),-!!groupname2) %>%
      rename(`sample name` = anno,
             `piRNA length` = pilen,
             `count of mismatch` = mismatch,
             `target length` = tag,
             `count of chimera` = count,
             `count of TE` = TEcount)
  }else{
    Vienna_stat_all <- bind_cols(refdata,Vienna_stat,Vienna_stat2,Vienna_stat3) %>%
      mutate(pilen = paste0(mean_pilen,"±",sd_pilen),
             mismatch = paste0(mean_mismatch,"±",sd_mismatch),
             tag = paste0(mean_tag,"±",sd_tag),
             dG = paste0(mean_dG,"±",sd_dG)) %>%
      select(-contains("mean"), -contains("sd"),-!!groupname2) %>%
      rename(`total reads` = comp,
             `end-to-end filtered reads` = comp2,
             `sample name` = anno,
             `piRNA length` = pilen,
             `count of mismatch` = mismatch,
             `target length` = tag,
             `count of chimera` = count,
             `count of TE` = TEcount)
  }
  Vienna_stat_all %>% gt() %>%
    gtsave(paste0(outdir,out,".pdf"))
  Vienna_stat_all %>% gt() %>%
    gtsave(paste0(outdir,out,".png"))
  Vienna_stat_all %>% gt() %>%
    gtsave(paste0(outdir,out,".html"))
  
  Vienna_stat_all %>% write_csv(paste0(outdir,out,".csv"))
  
  Vienna_stat_all2 <- Vienna_stat_all %>% as.matrix() %>% t() %>% as_tibble(rownames = NA)
  names(Vienna_stat_all2) <- Vienna_stat_all$`sample name`
  Vienna_stat_all2 <- Vienna_stat_all2  %>%
    rownames_to_column(var = "description") %>%
    filter(description != "sample name")
  Vienna_stat_all2 %>% write_csv(paste0(outdir,out,"2.csv"))
  
  Vienna_stat_all2 %>% gt() %>%
    gtsave(paste0(outdir,out,"2.pdf"))
  Vienna_stat_all2 %>% gt() %>%
    gtsave(paste0(outdir,out,"2.png"))
  Vienna_stat_all2 %>% gt() %>%
    gtsave(paste0(outdir,out,"2.html"))
}

Vienna_cluster_sle %>% filter(!(is.na(cluster))) %>% 
  make_table(refdata = tribble(~sample, ~group),groupname = cluster, out = "statscluster")

# fig2E
#######################################################################################################
# pairwise
#######################################################################################################

library(Biostrings)

test <- Vienna_cluster_sle %>%
  mutate(matchcount = str_count(piRNAmatch,"\\(")) %>%
  filter(pilen - matchcount < 5) %>%
  filter(str_detect(piRNAmatch, "\\.")) %>%
  select(piRNAseq,targetseq) %>%
  mutate(mate = map2(piRNAseq,targetseq, function(x,y) {
    pairwiseAlignment(pattern = x,
                      subject = reverse(y) %>%
                        chartr("CGAT", "GCTA",.),
                      type = "global-local",
                      gapOpening = -Inf) %>%
      compareStrings()
  }))


mismatch1 <- test %>% distinct() %>%
  mutate(miscount = str_count(mate, "\\?")) %>%
  filter(miscount == 1) %>%
  separate_rows(mate, sep = "") %>% filter(mate != "") %>%
  mutate(score = if_else(mate == "?",0,1)) %>%
  group_by(piRNAseq,targetseq) %>% mutate(site = row_number() %>% as.integer())
mismatch1 %>% filter(score ==0) %>%
  group_by(site) %>% count() %>% ungroup() %>% mutate(per = 100 * n/sum(n)) %>% 
  ggplot(aes(x = site, y = per))+
  geom_line(size = 1) +
  annotate("rect", xmin = 2,xmax = 8,ymin=0,ymax = Inf,alpha = 0.1, fill = "blue")+
  theme_minimal()+
  theme(panel.grid.minor.x =  element_blank(),panel.grid.minor.y =  element_blank())+
  xlab("piRNA base")+ylab("percentage of all bases (%)")+
  coord_cartesian(xlim = c(1, 30))+
  scale_x_continuous(breaks = c(1,5,10,15,20,25,30)) + scale_y_continuous(breaks = seq(0,7,1))+
  ggsave(paste0(outdir,"fig2E_mismatch1.png"),  width =3, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"fig2E_mismatch1.pdf"),  width =3, height =2.5)


mismatch2_con <- test %>% distinct() %>%
  mutate(miscount = str_count(mate, "\\?")) %>%
  filter(miscount == 2) %>% filter(str_detect(mate,"\\?\\?")) %>%
  separate_rows(mate, sep = "") %>% filter(mate != "") %>%
  mutate(score = if_else(mate == "?",0,1)) %>%
  group_by(piRNAseq,targetseq) %>% mutate(site = row_number() %>% as.integer())
mismatch2_con %>% filter(score ==0) %>% group_by(piRNAseq,targetseq) %>%
  summarise(site = mean(site)) %>% ungroup() %>% 
  group_by(site) %>% count() %>% ungroup() %>%mutate(per = 100 * n/sum(n)) %>% 
  ggplot(aes(x = site, y = per))+
  geom_line(size = 1) +
  annotate("rect", xmin = 2,xmax = 8,ymin=0,ymax = Inf,alpha = 0.1, fill = "blue")+
  theme_minimal()+
  theme(panel.grid.minor.x =  element_blank(),panel.grid.minor.y =  element_blank())+
  xlab("piRNA base")+ylab("percentage of all bases (%)")+
  coord_cartesian(xlim = c(1, 30))+
  scale_x_continuous(breaks = c(1,5,10,15,20,25,30)) + scale_y_continuous(breaks = seq(0,10,1))+
  ggsave(paste0(outdir,"fig2E_mismatch2_con.png"),  width =3, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"fig2E_mismatch2_con.pdf"),  width =3, height =2.5)


mismatch2 <- test %>% distinct() %>%
  mutate(miscount = str_count(mate, "\\?")) %>%
  filter(miscount == 2) %>% filter(!(str_detect(mate,"\\?\\?"))) %>%
  separate_rows(mate, sep = "") %>% filter(mate != "") %>%
  mutate(score = if_else(mate == "?",0,1)) %>%
  group_by(piRNAseq,targetseq) %>% mutate(site = row_number() %>% as.integer())
mismatch2 %>% filter(score ==0) %>%
  group_by(site) %>% count() %>% ungroup() %>% mutate(per = 100 * n/sum(n)) %>% 
  ggplot(aes(x = site, y = per))+
  geom_line(size = 1) +
  annotate("rect", xmin = 2,xmax = 8,ymin=0,ymax = Inf,alpha = 0.1, fill = "blue")+
  theme_minimal()+
  theme(panel.grid.minor.x =  element_blank(),panel.grid.minor.y =  element_blank())+
  xlab("piRNA base")+ylab("percentage of all bases (%)")+
  coord_cartesian(xlim = c(1, 30))+
  scale_x_continuous(breaks = c(1,5,10,15,20,25,30)) + scale_y_continuous(breaks = seq(0,7,1))+
  ggsave(paste0(outdir,"fig2E_mismatch2.png"),  width =3, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"fig2E_mismatch2.pdf"),  width =3, height =2.5)

mismatch3 <- test %>% distinct() %>%
  mutate(miscount = str_count(mate, "\\?")) %>%
  filter(miscount == 3) %>% filter(!(str_detect(mate,"\\?\\?\\?"))) %>%
  separate_rows(mate, sep = "") %>% filter(mate != "") %>%
  mutate(score = if_else(mate == "?",0,1)) %>%
  group_by(piRNAseq,targetseq) %>% mutate(site = row_number() %>% as.integer())
mismatch3 %>% filter(score ==0) %>%
  group_by(site) %>% count() %>% ungroup() %>% mutate(per = 100 * n/sum(n)) %>% 
  ggplot(aes(x = site, y = per))+
  geom_line(size = 1) +
  annotate("rect", xmin = 2,xmax = 8,ymin=0,ymax = Inf,alpha = 0.1, fill = "blue")+
  theme_minimal()+
  theme(panel.grid.minor.x =  element_blank(),panel.grid.minor.y =  element_blank())+
  xlab("piRNA base")+ylab("percentage of all bases (%)")+
  coord_cartesian(xlim = c(1, 30))+
  scale_x_continuous(breaks = c(1,5,10,15,20,25,30)) + scale_y_continuous(breaks = seq(0,7,1))+
  ggsave(paste0(outdir,"fig2E_mismatch3.png"),  width =3, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"fig2E_mismatch3.pdf"),  width =3, height =2.5)

#######################################################################################################
# Extended fig RNaseT1 G clustering
#######################################################################################################

Vienna_cluster_sle <- read_csv(paste0(outdir, "clash_sle_cluster3.csv"), col_types = cols(dG = col_double()))
bed_5base <- read_tsv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/","CLASH_target_5G.tsv"))

bed_5base2 <- bed_5base %>% 
  left_join(Vienna_cluster_sle %>% select(readname, sample,cluster), by = c("readname","sample")) %>% drop_na()

G5_logo_list <- bed_5base2 %>% select(cluster, site) %>%  group_nest(cluster) %>% 
  mutate(order = case_when(str_detect(cluster, "cluster1") ~ 1,
                           str_detect(cluster, "cluster2") ~ 2,
                           str_detect(cluster, "cluster3") ~ 3,
                           TRUE ~4)) %>% 
  arrange(order) %>% select(-order)
G5_logo_list2 <- map(1:length(G5_logo_list$cluster), function(x) pull(G5_logo_list$data[[x]]))
names(G5_logo_list2) <- G5_logo_list$cluster
ggseqlogo::ggseqlogo(G5_logo_list2,method = 'p',ncol = 1,seq_type = "rna")+
  annotate('rect',xmin = 5.5,xmax = 10.5,ymin = -0.05,ymax = 1,alpha = 0.5, fill = "gray")+
  annotate('segment',x = 5.5,xend = 10.5, y = -0.01,yend = -0.01,size = 1)+
  annotate("text",x = 8,y = -0.06,label = "target RNA")+
  scale_x_continuous(breaks=1:10, label = c("-5","-4","-3","-2","-1","1","2","3","4","5")) +
  ggsave(paste0(outdir,"CLASH_target_5G_cluster_probability.png"), width =3, height =7, dpi = 300)+
  ggsave(paste0(outdir,"CLASH_target_5G_cluster_probability.pdf"), width =3, height =7)


G3_logo_list <- Vienna_cluster_sle %>% 
  mutate(site = str_sub(targetseq,start = targetlen - 9 ,end=targetlen) %>%
           str_to_upper()%>% str_replace_all("T", "U")) %>%
  select(cluster, site) %>% drop_na() %>%  group_nest(cluster) %>% 
  mutate(order = case_when(str_detect(cluster, "cluster1") ~ 1,
                           str_detect(cluster, "cluster2") ~ 2,
                           str_detect(cluster, "cluster3") ~ 3,
                           TRUE ~4)) %>% 
  arrange(order) %>% select(-order)
G3_logo_list2 <- map(1:length(G3_logo_list$cluster), function(x) pull(G3_logo_list$data[[x]]))
names(G3_logo_list2) <- G3_logo_list$cluster
ggseqlogo::ggseqlogo(G3_logo_list2,method = 'p',ncol = 1,seq_type = "rna")+
  annotate('rect',xmin = 0.5,xmax = 10.5,ymin = -0.05,ymax = 1,alpha = 0.5, fill = "gray")+
  annotate('segment',x = 0.5,xend = 10.5, y = -0.01,yend = -0.01,size = 1)+
  annotate("text",x = 5,y = -0.06,label = "target RNA")+
  scale_x_continuous(breaks=1:10, label = c("-10","-9","-8","-7","-6","-5","-4","-3","-2","-1")) +
  ggsave(paste0(outdir,"CLASH_target_last5G_cluster_probability.png"), width =3, height =7, dpi = 300)+
  ggsave(paste0(outdir,"CLASH_target_last5G_cluster_probability.pdf"), width =3, height =7)

# piRNA last G

G3_logo_list3 <- Vienna_cluster_sle %>% 
  mutate(site = str_sub(piRNAseq,start = pilen - 9 ,end=pilen) %>%
           str_to_upper()%>% str_replace_all("T", "U")) %>%
  select(cluster, site) %>% drop_na() %>%  group_nest(cluster) %>% 
  mutate(order = case_when(str_detect(cluster, "cluster1") ~ 1,
                           str_detect(cluster, "cluster2") ~ 2,
                           str_detect(cluster, "cluster3") ~ 3,
                           TRUE ~4)) %>% 
  arrange(order) %>% select(-order)
G3_logo_list4 <- map(1:length(G3_logo_list3$cluster), function(x) pull(G3_logo_list3$data[[x]]))
names(G3_logo_list4) <- G3_logo_list3$cluster
ggseqlogo::ggseqlogo(G3_logo_list4,method = 'p',ncol = 1,seq_type = "rna")+
  annotate('rect',xmin = 0.5,xmax = 10.5,ymin = -0.05,ymax = 1,alpha = 0.5, fill = "gray")+
  annotate('segment',x = 0.5,xend = 10.5, y = -0.01,yend = -0.01,size = 1)+
  annotate("text",x = 5,y = -0.06,label = "piRNA")+
  scale_x_continuous(breaks=1:10, label = c("-10","-9","-8","-7","-6","-5","-4","-3","-2","-1")) +
  ggsave(paste0(outdir,"CLASH_target_pi_last5G_cluster_probability.png"), width =3, height =7, dpi = 300)+
  ggsave(paste0(outdir,"CLASH_target_pi_last5G_cluster_probability.pdf"), width =3, height =7)

