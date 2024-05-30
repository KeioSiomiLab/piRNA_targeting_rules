

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/CLASH/vienna/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gplots)))
suppressMessages(suppressWarnings(require(RColorBrewer)))
suppressMessages(suppressWarnings(require(gtools)))
suppressMessages(suppressWarnings(require(gridExtra)))
suppressMessages(suppressWarnings(require(stringi)))
suppressMessages(suppressWarnings(require(gt)))
suppressMessages(suppressWarnings(require(stringi)))
ppi <- 600
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

outdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/fig1/"
CLASHdir <- "/media/pericles/CLASH/"
CLASH_shufdir <- "/media/pericles/CLASH/shuffle/"


#######################################################################################################
# table making
#######################################################################################################
piRNA_class <- read_tsv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/piRNA/", "piRNA_classification.tsv")) %>% 
  rename(piRNA = name, piRNAtype = type3)

Vienna_sle <-read_csv(paste0(CLASHdir, "vienna1fig/clash_colapse.csv"),
                      col_types = cols(dG = col_double())) %>% filter(group!= "siHP1a") %>% 
  select(-piRNAtype) %>% left_join(piRNA_class, by = "piRNA")
Vienna_sle %>% write_csv(paste0(outdir,"CLASH_sup_file1.csv"))

Vienna_select <- read_csv(paste0(CLASHdir, "vienna1fig/clash_selection.csv"), col_types = cols(dG = col_double())) %>% 
  filter(group!= "siHP1a") %>% select(-piRNAtype) %>% left_join(piRNA_class, by = "piRNA")
Vienna_select %>% write_csv(paste0(outdir,"CLASH_sup_file2.csv"))

Vienna_shuffle <- read_tsv(paste0(CLASH_shufdir, "Vienna_sle_vienna_shuffle.tsv"), col_types = cols(dG = col_double())) %>% 
  filter(!(str_detect(sample,"HP1")))
Vienna_shuffle %>% write_csv(paste0(outdir,"CLASH_sup_file3.csv"))

Vienna_struct <- read_csv(paste0(CLASHdir, "vienna1fig/clash_selection_struct.csv"), col_types = cols(dG = col_double()))  %>% 
  filter(group!= "siHP1a") %>% select(-piRNAtype) %>% left_join(piRNA_class, by = "piRNA")
Vienna_struct %>% write_csv(paste0(outdir,"CLASH_sup_file4.csv"))

# tag-pi shuffle
CLASH_shufdir_tag <- "/media/pericles/CLASH/shuffle_tag/"
Vienna_shuffle_tag <- read_tsv(paste0(CLASH_shufdir_tag, "Vienna_sle_vienna_shuffle_tag.tsv"), col_types = cols(dG = col_double())) %>% 
  filter(!(str_detect(sample,"HP1")))
Vienna_shuffle_tag %>% write_csv(paste0(outdir,"CLASH_sup_file3_tag.csv"))

# wobble pair annotation
Vienna_sle  <- read_csv(paste0(outdir,"CLASH_sup_file1.csv"), col_types = cols(dG = col_double()))

wobble_pi <- Vienna_sle  %>% select(readname, sample, piRNAseq, piRNAmatch) %>% 
  mutate(piRNAmatch2 = piRNAmatch %>% str_replace_all("\\.","0") %>% str_replace_all("\\(","1") %>% str_split(pattern = ""),
         piRNAseq2 = piRNAseq %>% str_split(pattern = ""), index2 = row_number()) %>% 
  mutate(base = map(piRNAmatch2, function(x) {seq(1, length(x), by = 1)})) %>% 
  select(sample, readname, index2, base, piRNAmatch2, piRNAseq2) %>% 
  unnest_legacy() %>% filter(piRNAmatch2 == "1") %>% 
  group_by(sample, readname, index2) %>% mutate(base2 = row_number()) %>% ungroup() %>% select(-base, -piRNAmatch2)

wobble_tag <- Vienna_sle %>% select(readname, sample, targetseq, targetmatch) %>% 
  mutate(targetmatch2 = targetmatch %>% stringi::stri_reverse() %>% str_replace_all("\\.","0") %>% str_replace_all("\\)","1") %>% str_split(pattern = ""),
         targetseq2 = targetseq %>% stringi::stri_reverse() %>% str_split(pattern = ""), index2 = row_number()) %>% 
  mutate(base = map(targetmatch2, function(x) {seq(1, length(x), by = 1)})) %>% 
  select(sample, readname, index2, base, targetmatch2, targetseq2) %>% 
  unnest_legacy() %>% filter(targetmatch2 == "1") %>% 
  group_by(sample, readname, index2) %>% mutate(base2 = row_number())%>% ungroup() %>% select(-base, -targetmatch2)

wobble_all <- full_join(wobble_pi, wobble_tag, by = c("sample", "readname", "index2","base2")) %>% 
  mutate(wobble = if_else((piRNAseq2 == "G" & targetseq2 == "T") | (piRNAseq2 == "T" & targetseq2 == "G"), "wobble","pair" )) %>% 
  filter(wobble == "wobble") %>% select(readname, sample) %>% distinct()
wobble_all %>% write_csv(paste0(outdir,"CLASH_wobble_annotation.csv"))
#######################################################################################################
# fig1B&C&D&E
#######################################################################################################

# fig1B&C&D
Vienna_shuffle <- read_csv(paste0(outdir, "CLASH_sup_file3.csv"), col_types = cols(dG = col_double())) %>%
  mutate(dG = if_else(!(str_detect(piRNAmatch, "\\(")), Inf, dG)) %>%
  select(sample, readname,dG,piRNAmatch,targetmatch) %>%
  rename(dG_shuf = dG ,
         piRNAmatch_shuf = piRNAmatch,
         targetmatch_shuf = targetmatch) %>%
  right_join(Vienna_sle, by = c("sample","readname"))

Vienna_shuffle %>% filter(group!= "siHP1a") %>% select(dG_shuf,dG,sample, readname,group) %>%
  rename(shuffle = dG_shuf, chimera = dG) %>%
  pivot_longer(cols = shuffle:chimera) %>%
  mutate(value = if_else(is.infinite(value),0,value),
         binneddG = round(value)) %>% 
  group_by(binneddG, name) %>% summarise(count =n()) %>% 
  ggplot(aes(x = binneddG, color = name, y = count)) +
  geom_line(size=1) + geom_area(aes(fill = name),alpha = 0.1, position = "identity") +
  guides(color = "none", fill = "none")+
  labs(x = "dG (kcal/mol)", y = "Number of chimeric reads") +
  theme_minimal()+ 
  ggsave(paste0(outdir,"fig1D_dG_collapse.png"), width =3, height =3, dpi = 300) +
  ggsave(paste0(outdir,"fig1D_dG_collapse.pdf"), width =3, height =3)
  
Vienna_shuffle2 <- Vienna_shuffle %>% select(sample,readname,dG,piRNAmatch) %>%
  mutate(type = "chimera") %>%
  bind_rows(Vienna_shuffle %>% select(sample,readname,dG_shuf,piRNAmatch_shuf) %>%
              rename(dG =dG_shuf,piRNAmatch=piRNAmatch_shuf) %>%
              mutate(type = "shuffle")) %>%
  rename(match = piRNAmatch)

Vienna_shuffle2 %>% mutate(mismatch = str_count(match, "\\."),
                           group = if_else(str_detect(sample, "HP1"),"siHP1a", "WT"),
                           type = if_else(type == "chimera", "CLASH-chimeras", "shuffled-chimeras")) %>%
  filter(group!= "siHP1a") %>%
  group_by(type, mismatch) %>%
  summarise(count = n()) %>% ungroup() %>%
  ggplot(aes(x = mismatch, y = count, color = type))+
  geom_line(size=1)+ geom_area(aes(fill = type),alpha = 0.1, position = "identity", color = NA) +
  labs(color = "", y = "Number of chimeric reads", x="mismatch number") +
  theme_minimal() + guides(color=guide_legend(override.aes=list(fill=NA)), fill = "none")+
  theme(legend.position = c(0.7,0.7),legend.key=element_blank())+
  ggsave(paste0(outdir,"fig1D_mismatch_collapse.png"), width =3, height =3, dpi = 300) +
  ggsave(paste0(outdir,"fig1D_mismatch_collapse.pdf"), width =3, height =3)

cal_heat <- function(data) {
  heatmap <- data %>%
    mutate(heat = match %>% str_replace_all("\\.","0-") %>% str_replace_all("\\(","1-") %>% str_replace_all("\\)","1-"),
           dG = if_else(is.infinite(dG), 0,dG)) %>%
    mutate(score = map(heat, function(x) {x %>% str_sub(end = -2) %>% str_split(pattern = "-")}),
           index2 = row_number()) %>%
    unnest(cols = score) %>%
    mutate(base = map(score, function(x) {seq(1, length(x), by = 1)})) %>%
    select(sample,index2, dG,score, base,readname,type) %>%
    unnest_legacy() %>% drop_na() %>%
    mutate(score_G = as.integer(score) * dG) %>%
    unite(sample,readname, sep="&",col="gene_name")
  return(heatmap)
}
heatmap_cal <- function(data, length = 30,out){
  data %>% filter(base <= length) %>%
    ggplot() +
    geom_raster(aes(y =fct_reorder(gene_name,dG) %>% fct_rev(), fill = score_G, x=base))+
    facet_wrap(~type, strip.position ="top", scales = "free_y", ncol=2)+
    scale_x_continuous(breaks = c(1,10,20))+theme_minimal()+
    scale_fill_gradient2(low = "black",mid = "black", high = "white", midpoint = -50)+
    labs(fill = "\u0394G (kcal/mol)", x = "bases (nt)", y= "")+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          strip.background=element_blank())+
    ggsave(paste0(outdir,"fig1B_collapse.png"), width =3.5, height =6, dpi = 600)
}

Vienna_shuffle3 <- Vienna_shuffle2 %>% 
  mutate(mismatch = str_count(match, "\\."),
         group = if_else(str_detect(sample, "HP1"),"siHP1a", "WT")) %>%
  filter(group!= "siHP1a") %>% 
  cal_heat()
Vienna_shuffle3 %>% mutate(type = if_else(type=="chimera","CLASH-chimeras","shuffled-chimeras")) %>% 
  heatmap_cal(length = 30,out = "piRNA")

Vienna_shuffle4 <- Vienna_shuffle2 %>% mutate(mismatch = str_count(match, "\\."),
                                              group = if_else(str_detect(sample, "HP1"),"siHP1a", "WT")) %>%
  filter(group!= "siHP1a") %>% mutate(pilen=str_length(match)) %>%
  mutate(heat = match %>% str_replace_all("\\.","0-") %>% str_replace_all("\\(","1-") %>% str_replace_all("\\)","1-"),
         dG = if_else(is.infinite(dG), 0,dG)) %>%
  mutate(score = map(heat, function(x) {x %>% str_sub(end = -2) %>% str_split(pattern = "-")}),
         index2 = row_number()) %>%
  unnest(cols = score) %>%
  mutate(base = map(score, function(x) {seq(1, length(x), by = 1)})) %>%
  select(sample,index2, dG,score, base,readname,type,  pilen) %>%
  unnest_legacy()

Vienna_shuffle4 %>%
  filter(pilen %in% c(24,25,26,27,28)) %>%
  mutate(readcount = readname %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer(),
         score = score %>% as.double(),
         pilen = paste0(pilen, "nt")) %>%
  group_by(type,pilen, base) %>%
  summarise(base_pair = sum(score * readcount)/ sum(readcount)) %>% ungroup()%>%
  ggplot(aes(x = base, y = 100 *base_pair, color =pilen)) +
  geom_line(aes(linetype = type)) +
  annotate("rect", xmin = 2,xmax = 8,ymin=0,ymax = 100,alpha = 0.1, fill = "blue")+
  theme_minimal()+
  coord_cartesian(xlim = c(1, 28), ylim = c(0,100)) +
  scale_x_continuous(breaks = c(1,5,10,15,20,25,28)) +theme(panel.grid.minor.x = element_blank())+
  labs(color = "", x = "bases (nt)", y= "base pairing (%)", linetype ="") +
  ggsave(paste0(outdir,"fig1C.png"),  width =4, height =3, dpi = 300)+
  ggsave(paste0(outdir,"fig1C.pdf"),  width =4, height =3)

#######################################################################################################
# fig1E
#######################################################################################################

referenceRNAorder <- c("TE","mRNA","ncRNA","rRNA","snoRNA","snRNA","tRNA","pseudogene")
# fig1E
Vienna_sle %>% mutate(targettype = if_else(targettype %in% c("mRNA(TE)","ensemblTE","rmblastTE","pbsvTE"), "TE",targettype)) %>%
  filter(group!= "siHP1a") %>% 
  mutate(readcount = readname %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer()) %>% 
  group_by(targettype) %>%
  summarise(n = sum(readcount)) %>% ungroup() %>%  mutate(prop = 100 *n/sum(n)) %>% 
  mutate(ypos = cumsum(prop)-0.5*prop) %>%
  mutate(targettype = targettype %>% factor(levels = referenceRNAorder %>% fct_rev())) %>% 
  ggplot(aes(x="", y=prop, fill=targettype)) +
  geom_bar(stat="identity", width=1, color="white", linewidth=0.1) +
  coord_polar("y", start=0, direction = 1) +
  theme_void() +
  guides(fill= guide_legend(title = NULL))+
  geom_text(aes(label = if_else(prop > 3,paste0(round(prop,1), "%"),"")), position = position_stack(vjust = 0.5))+
  ggsave(paste0(outdir,"fig1E_pie.png"), width =4, height =3, dpi = 300)+
  ggsave(paste0(outdir,"fig1E_pie.pdf"), width =4, height =3)

RNAtable_CLASH <- Vienna_sle %>% mutate(targettype = if_else(targettype %in% c("mRNA(TE)","ensemblTE","rmblastTE","pbsvTE"), "TE",targettype)) %>%
  filter(group!= "siHP1a") %>% 
  mutate(readcount = readname %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer()) %>% 
  group_by(targettype) %>%
  summarise(n = sum(readcount)) %>% ungroup() %>%  mutate(prop = 100 *n/sum(n)) %>% 
  mutate(ypos = cumsum(prop)-0.5*prop)
RNAtable_CLASH %>% 
  write_csv(paste0(outdir,"CLASH_RNAtype_table.csv"))

#######################################################################################################
# fig1F
#######################################################################################################
comparedir <- "/media/pericles/CLASH/compare/"

files3 <- dir(comparedir , pattern = "_struct.csv$", full.names = TRUE)
files4 <- files3 %>%
  str_replace_all("_struct.csv$", "") %>%
  str_replace(paste0(comparedir, "/"), "")
df2 <- vector("list", length(files3))
for (i in seq_along(files3)) {
  df2[[i]] <- read_tsv(files3[[i]], col_types = cols(dG = col_double())) %>%
    select(sample,readname,dG,piRNAseq,targetseq,piRNAmatch,targetmatch,ct)
  df2[[i]][["dataset"]] <- files4[[i]]
}
compare_vienna <- bind_rows(df2) %>%
  mutate(struct = if_else(str_detect(ct,"interior"),"interior",
                          if_else(str_detect(ct,"bulge"),"bulge",
                                  if_else(str_detect(ct,"mismatch"),"mismatch",ct)))) %>%
  mutate(struct=if_else(struct=="perfect" & str_detect(piRNAmatch,"\\."),"mate",struct)) %>% 
  filter(dataset!="hiCLIP") %>% 
  mutate(dataset = case_when(dataset == "AGO" ~ "hAGO1",
                             dataset == "Aub" ~ "dAub",
                             TRUE ~ dataset)) %>% 
  filter(!(str_detect(sample, "HP1_"))) 
sample_order1 <- c("Piwi","hAGO1", "dAub", "PRG1", "SMEDWI3")
sample_order2 <- c("perfect","mate","mismatch","bulge","interior")

compare_vienna %>% mutate(dataset = dataset %>% factor(levels = sample_order1)) %>% 
  group_by(dataset) %>% count() %>% ungroup() %>% 
  ggplot(aes(x = n,y = dataset)) +
  geom_col(fill = "orange")+ geom_text(aes(x = 120000,label = paste0("n = ",n)), size = 2.5)+
  labs(x = "detected chimera count")+  coord_cartesian(xlim = c(0, 130000)) +
  theme_minimal()+
  theme(axis.title.y = element_blank()) +
  ggsave(paste0(outdir,"dataset_count.png"), width =4, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"dataset_count.pdf"), width =4, height =2.5)

compare_vienna %>% mutate(dataset = dataset %>% factor(levels = sample_order1)) %>% filter(struct !="else") %>% 
  mutate(dG = dG %/% 1) %>%
  group_by(dataset, dG) %>% summarise(count = n()) %>% ungroup() %>%
  group_by(dataset) %>% mutate(per = 100 * count/ sum(count)) %>% ungroup() %>% 
  ggplot(aes(x = dG, y = per, fill = dataset, color = dataset)) +
  geom_line()+geom_area(alpha = 0.1, position = "identity", color= NA)+ #geom_point()+
  theme_minimal()+ 
  theme(legend.position = "bottom", legend.title = element_blank(),legend.key=element_blank()) +
  guides(color=guide_legend(override.aes=list(fill=NA),nrow = 2),fill = "none")+
  xlab("dG (kcal/mol)") + ylab("percentage of all pair (%)") + 
  ggsave(paste0(outdir,"Ext_other_exp_compare_dG3.png"), width =4, height =3.3, dpi = 300)+
  ggsave(paste0(outdir,"Ext_other_exp_compare_dG3.pdf"), width =4, height =3.3)

compare_vienna %>% mutate(mismatch = str_count(piRNAmatch,"\\.")) %>% 
  group_by(dataset, mismatch) %>% summarise(count = n()) %>% ungroup() %>%
  group_by(dataset) %>% mutate(per = 100 * count/ sum(count)) %>% ungroup() %>% 
  mutate(dataset = dataset %>% factor(levels = sample_order1)) %>% 
  ggplot(aes(x = mismatch, y = per, color = dataset, fill = dataset))+ 
  geom_line()+geom_area(alpha = 0.1, position = "identity", color= NA)+ #geom_point()+
  xlab("mismatch count") + ylab("percentage of all pair (%)") +theme_minimal() +
  theme(legend.position = "bottom",legend.title = element_blank(),legend.key=element_blank())+
  guides(color=guide_legend(override.aes=list(fill=NA),nrow = 2),fill = "none")+
  ggsave(paste0(outdir,"Ext_other_exp_compare_mismatch.png"), width =4, height =3.3, dpi = 300)+
  ggsave(paste0(outdir,"Ext_other_exp_compare_mismatch.pdf"), width =4, height =3.3)

total <- compare_vienna %>% group_by(dataset) %>%
  summarise(total = n()) %>% ungroup()
compare_vienna2 <- compare_vienna %>% group_by(dataset,struct) %>%
  summarise(count = n()) %>% ungroup() %>%
  left_join(total,by="dataset") %>%
  mutate(per = 100*count/total)

compare_vienna2 %>% filter(struct !="else") %>% 
bind_rows(tibble(dataset ="dAub",
                   struct = "perfect",
                   count = 0,
                   total = 31852,
                   per =0)) %>% 
  mutate(dataset = dataset %>% factor(levels = sample_order1),
         struct = struct %>% factor(levels = sample_order2)) %>% 
  ggplot(aes(x = dataset,fill=per,y = struct %>% fct_rev()))+
  geom_tile()+
  geom_text(aes(label = per %>% round(digits = 2)),
            colour = "black",size = 3)+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 30, hjust=0.8, vjust=1))+
  labs(x = "", y = "", fill = "% chimera")+
  scale_fill_gradient(low="white", high="dodgerblue3") +
  ggsave(paste0(outdir,"fig1F.png"), width =3.5, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"fig1F.pdf"), width =4, height =2.5)



#######################################################################################################
# Supplementary figure
#######################################################################################################
# Sup fig1B
library(GGally)
library(ggthemes)

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


Vienna_select <- read_csv(paste0(outdir,"CLASH_sup_file2.csv"), col_types = cols(dG = col_double()))
target_scatter <- Vienna_select %>% filter(targettype %in% c("mRNA","ncRNA","ensemblTE")) %>% 
  select(sample, targetname,readname, targettype) %>%
  group_by(sample, targetname,targettype) %>% summarise(count = n()) %>% ungroup()
target_scatter2 <- target_scatter %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count = 0))

png(paste0(outdir,"sup_fig1B.png"), width = 3*ppi, height = 3*ppi, res = ppi)
pm <- ggpairs(data = target_scatter2 %>% select(-contains("target")) %>%
                mutate_at(vars(contains("rep")),list(~log10(.+1))),
              xlab = "log10(target count + 1) ",ylab = "log10(target count + 1)",
              lower = list(continuous = wrap("points", alpha = 0.3, size = 0.1)),
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

#######################################################################################################
# Sup fig1D
#######################################################################################################

Vienna_sle %>% group_by(pilen) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = pilen, y = count)) + geom_col(fill = "orange")+
  xlab("piRNA length") + ylab("Number of chimeric reads") +
  geom_text(aes(label = count), vjust = -0.5, hjust=0.5, colour = "black",size = 3)+
  theme_minimal() +
  scale_x_continuous(breaks = c(23,24,25,26,27,28,29,30))+
  ggsave(paste0(outdir,"sup_fig1D_colapse.png"), width =4, height =3.5, dpi = 300)+
  ggsave(paste0(outdir,"sup_fig1D_colapse.pdf"), width =4, height =3.5)

#######################################################################################################
# Extended fig RNaseT1 G
#######################################################################################################
Vienna_sle <-read_csv(paste0(outdir, "CLASH_sup_file1.csv"), col_types = cols(dG = col_double()))
Vienna_select <- read_csv(paste0(outdir,"CLASH_sup_file2.csv"), col_types = cols(dG = col_double()))

targetRNAlevel <- c("rRNA", "tRNA", "snRNA", "snoRNA", "ensemblTE","pbsvTE","rmblastTE", "mRNA(TE)","mRNA",
                    "pseudogene", "ncRNA", "pre_miRNA", "miRNA", NA)

Vienna_RNaseT1 <- Vienna_sle %>% 
  select(sample, readname,dG, piRNA,piRNAtype,pilen, piRNAseq, targettype) %>% 
  left_join(Vienna_select %>% select(readname, sample,target, bedstart,bedend), by = c("readname","sample")) %>% 
  group_by(sample, readname) %>%
  slice_max(n=1, target, with_ties = FALSE) %>% ungroup()
  
gene_seq <- read_tsv(paste0("/media/pericles/CLASH/","database/dm6_gene.tsv"),
                     col_names = c("target","seq")) %>%
  separate(target,into = c("target","targetname","targettype"), sep="_") %>%
  mutate(target = target %>% str_remove("\\-")) %>%
  select(target, seq)

bed_5base <- Vienna_RNaseT1 %>% 
  mutate(new_base = bedstart -5+1,
         new_base2 = bedstart+5) %>%
  filter(new_base > 0) %>%
  left_join(gene_seq,by = "target") %>%
  mutate(site = str_sub(seq,start = new_base,end=new_base2) %>%
           str_to_upper()%>% str_replace_all("T", "U")) %>%
  filter(!(str_detect(site,"N"))) %>% 
  select(sample, readname,site)

ggseqlogo::ggseqlogo(bed_5base$site,method = 'p',seq_type = "rna")+
  annotate('rect',xmin = 5.5,xmax = 10.5,ymin = -0.05,ymax = 1,alpha = 0.5, fill = "gray")+
  annotate('segment',x = 5.5,xend = 10.5, y = -0.01,yend = -0.01,size = 1)+
  annotate("text",x = 8,y = -0.06,label = "target RNA")+
  scale_x_continuous(breaks=1:10, label = c("-5","-4","-3","-2","-1","1","2","3","4","5")) +
  ggsave(paste0(outdir,"CLASH_target_5G_probability.png"), width =3, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"CLASH_target_5G_probability.pdf"), width =3, height =2.5)

bed_5base %>% write_tsv(paste0(outdir,"CLASH_target_5G.tsv"))

Vienna_RNaseT1_2 <- Vienna_sle %>% 
  mutate(site = str_sub(targetseq,start = targetlen - 9 ,end=targetlen) %>%
           str_to_upper()%>% str_replace_all("T", "U")) %>%
  select(sample, readname,site)

ggseqlogo::ggseqlogo(Vienna_RNaseT1_2$site,method = 'p',seq_type = "rna")+
  annotate('rect',xmin = 0.5,xmax = 10.5,ymin = -0.05,ymax = 1,alpha = 0.5, fill = "gray")+
  annotate('segment',x = 0.5,xend = 10.5, y = -0.01,yend = -0.01,size = 1)+
  annotate("text",x = 5,y = -0.06,label = "target RNA")+
  scale_x_continuous(breaks=1:10, label = c("-10","-9","-8","-7","-6","-5","-4","-3","-2","-1")) +
  ggsave(paste0(outdir,"CLASH_target_last5G_probability.png"), width =3, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"CLASH_target_last5G_probability.pdf"), width =3, height =2.5)

# piRNA last sequence
Vienna_RNaseT1_3 <- Vienna_sle %>% 
  mutate(site = str_sub(piRNAseq,start = pilen - 9 ,end=pilen) %>%
           str_to_upper()%>% str_replace_all("T", "U")) %>%
  select(sample, readname,site)

ggseqlogo::ggseqlogo(Vienna_RNaseT1_3$site,method = 'p',seq_type = "rna")+
  annotate('rect',xmin = 0.5,xmax = 10.5,ymin = -0.05,ymax = 1,alpha = 0.5, fill = "gray")+
  annotate('segment',x = 0.5,xend = 10.5, y = -0.01,yend = -0.01,size = 1)+
  annotate("text",x = 5,y = -0.06,label = "piRNA")+
  scale_x_continuous(breaks=1:10, label = c("-10","-9","-8","-7","-6","-5","-4","-3","-2","-1")) +
  ggsave(paste0(outdir,"CLASH_target_pi_last5G_probability.png"), width =3, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"CLASH_target_pi_last5G_probability.pdf"), width =3, height =2.5)

#######################################################################################################
# CLASH stats
#######################################################################################################


read_count1 <- read_tsv("/media/pericles/CLASH/fasta_tsv/comp_stats.tsv")
samplename <- c("rep1", "rep2","rep3","rep4","rep5","HP1_rep1","HP1_rep2","HP1_rep3")
read_count2 <- read_count1 %>%
  filter(str_detect(file, "rep")) %>%
  filter(str_detect(file,"comp|comp2")) %>%
  filter(!(str_detect(file,"comp_tf"))) %>%
  filter(!(str_detect(file,"test"))) %>%
  separate(file, into = c("dis", "name"),sep = "\\/") %>%
  mutate(name = name %>% str_remove(".fastq|.fasta"),
         group = if_else(str_detect(name,"HP1"), "siHP1a", "WT"),
         name2= name %>% str_remove("HP1_")) %>%
  separate(name2,into = c("sample","procedure"),sep="_") %>%
  mutate(anno=paste0(group, "_", sample)) %>% 
  filter(group == "WT")

read_count3 <- read_count2 %>%
  mutate(name = name %>% str_remove("_comp$") %>% str_remove("_comp2$")) %>%
  select(name,procedure,anno,num_seqs) %>%
  pivot_wider(names_from = procedure,values_from = num_seqs) %>%
  filter(name !="rep6") %>%
  select(-name)

read_count3_KD <- read_count2 %>%
  mutate(name = name %>% str_remove("_comp$") %>% str_remove("_comp2$")) %>%
  select(group,name,procedure,anno,num_seqs) %>%
  pivot_wider(names_from = procedure,values_from = num_seqs) %>%
  filter(name !="rep6") %>% 
  group_by(group) %>% summarise(comp = sum(comp),comp2 = sum(comp2)) %>%
  ungroup() %>% rename(anno = group)


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

Vienna_sle %>% filter(group == "WT") %>% 
  make_table(refdata = read_count3, groupname = sample, out = "stats")

Vienna_sle %>% filter(group == "WT") %>% 
  make_table(refdata = read_count3_KD,groupname = group, out = "statsKD")


#######################################################################################################
# pi-tag & tag-pi
#######################################################################################################
# Ext data

Vienna_sle_pitag <- read_csv(paste0(CLASHdir, "vienna1fig/clash_colapse_all.csv"), 
                             col_types = cols(dG = col_double())) %>% filter(group!= "siHP1a") %>% 
  mutate(direction = if_else(direction == "pi_tag", "pi-tag","tag-pi")) %>% 
  mutate(readcount = readname %>% str_extract("_\\d+") %>% str_remove("_") %>% as.integer())

Vienna_sle_pitag %>% group_by(direction) %>% count() %>%ungroup() %>% 
  ggplot(aes(x = n,y = direction)) +
  geom_col(fill = "orange")+ geom_text(aes(x = 20000,label = paste0("n = ",n)), size = 2.5)+
  labs(x = "detected chimera count")+  coord_cartesian(xlim = c(0, 22000)) +
  theme_minimal()+ theme(axis.title.y = element_blank()) +
  ggsave(paste0(outdir,"Ext_direction_count.png"), width =4, height =1.5, dpi = 300)+
  ggsave(paste0(outdir,"Ext_direction_count.pdf"), width =4, height =1.5)

Vienna_shuffle <-  read_csv(paste0(outdir,"CLASH_sup_file3.csv"),col_types = cols(dG = col_double())) %>% 
  mutate(direction = "pi-tag") %>% select(sample, readname, direction,dG) %>% mutate(group2 = "shuffle")
Vienna_shuffle_tag  <-  read_csv(paste0(outdir,"CLASH_sup_file3_tag.csv"),col_types = cols(dG = col_double())) %>% 
  mutate(direction = "tag-pi") %>% select(sample, readname, direction,dG) %>% mutate(group2 = "shuffle")

Vienna_sle_pitag %>% select(sample, readname, direction,dG) %>% mutate(group2 = "chimera") %>% bind_rows(Vienna_shuffle ,Vienna_shuffle_tag) %>% 
  mutate(group2 = if_else(group2 == "chimera", "CLASH-chimeras", "shuffled-chimeras")) %>% 
  ggplot(aes(x = dG, color = direction)) + 
  geom_freqpoly(aes(linetype = group2),binwidth=1) +
  theme_minimal()+ guides(color=guide_legend(nrow = 2,byrow = TRUE), linetype=guide_legend(nrow = 2,byrow = TRUE)) +
  theme(legend.position = "bottom", legend.title = element_blank(), )+
  labs(x = "dG (kcal/mol)", y = "Number of chimeric reads") +
  ggsave(paste0(outdir,"Ext_direction_distri_poly.png"), width =3.5, height =3.8, dpi = 300)+
  ggsave(paste0(outdir,"Ext_direction_distri_poly.pdf"), width =3.5, height =3.8)
#
referenceRNAorder <- c("TE","mRNA","ncRNA","rRNA","snoRNA","snRNA","tRNA","pseudogene")
Vienna_sle_pitag %>% mutate(targettype = if_else(targettype %in% c("mRNA(TE)","ensemblTE","rmblastTE","pbsvTE"), "TE",targettype)) %>%
  group_by(direction,group) %>%
  count(targettype) %>% mutate(per = 100 *n/sum(n)) %>% ungroup() %>%
  ggplot(aes(x = direction, y = per, fill = fct_reorder(targettype, per) %>% fct_rev())) +
  geom_col(color = "black")+
  geom_text(aes(label = if_else(round(per,1)>5, round(per,1) %>% as.character(), "")),
            colour = "black",size = 3, position = position_stack(vjust = 0.5))+
  theme_minimal() + coord_flip() +
  labs(fill = "RNA type", y = "Fraction (%)") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank()) +
  ggsave(paste0(outdir,"Ext_direction_target_RNAtype_all.png"), width =5, height =2.5, dpi = 300)+
  ggsave(paste0(outdir,"Ext_direction_target_RNAtype_all.pdf"), width =5, height =2.5)


#######################################################################################################
# supplementary table (compare dG)
#######################################################################################################
# Show the selecting strategy

Vienna1 <- read_csv(paste0("/media/pericles/CLASH/", "vienna1fig/clash.csv"), col_types = cols(dG = col_double())) %>% 
  select(sample, readname,dG, piRNA,piRNAtype,pilen, piRNAseq, targetlen,targetseq,piRNAmatch, targetmatch,group,direction) %>% 
  filter(direction=="pi_tag") %>%filter(group!= "siHP1a") %>% distinct()

Vienna1_double <- Vienna1 %>% 
  group_by(readname, sample) %>% mutate(total = n()) %>% ungroup() %>% filter(total!= 1)

Vienna1_test <- Vienna1 %>% 
  group_by(readname, sample) %>% mutate(total = n()) %>% ungroup() %>% filter(total== 5)

#increase by adding
# rep2, 1385815_1
Vienna1_double %>% filter(sample == "rep2" & readname == "1385815_1") %>% 
  select(sample, readname,dG, piRNA,pilen, piRNAseq, targetlen,targetseq,piRNAmatch, targetmatch) %>% 
  write_csv(paste0(outdir,"CLASH_Ext_decide_lig_case1_2.csv"))

# not change
# rep2, 1382882_1
Vienna1_double %>% filter(sample == "rep2" & readname == "1382882_1")  %>% 
  select(sample, readname,dG, piRNA,pilen, piRNAseq, targetlen,targetseq,piRNAmatch, targetmatch) %>% 
  write_csv(paste0(outdir,"CLASH_Ext_decide_lig_case2_2.csv"))

# mix
# rep1, 104986_2
Vienna1_double %>% filter(sample == "rep1" & readname == "104986_2") %>% 
  select(sample, readname,dG, piRNA,pilen, piRNAseq, targetlen,targetseq,piRNAmatch, targetmatch) %>% 
  write_csv(paste0(outdir,"CLASH_Ext_decide_lig_case3_1.csv"))


#######################################################################################################
# first unpaired
#######################################################################################################

tes <- Vienna_sle %>% filter(str_detect(piRNAmatch,"^\\.\\.\\(")) 
  
# detect 2mismatch and 1 unpaired
# rep1, 1473837_1
Vienna_sle %>% filter(sample == "rep1" & readname == "1473837_1") %>% 
  select(sample, readname,dG, piRNA,pilen, piRNAseq, targetlen,targetseq,piRNAmatch, targetmatch) %>% 
  write_csv(paste0(outdir,"CLASH_Ext_first_unpaired_1.csv"))
# rep1, 1492859_1
Vienna_sle %>% filter(sample == "rep1" & readname == "1492859_1") %>% 
  select(sample, readname,dG, piRNA,pilen, piRNAseq, targetlen,targetseq,piRNAmatch, targetmatch) %>% 
  write_csv(paste0(outdir,"CLASH_Ext_first_unpaired_2.csv"))

# detect 3mismatch and 1 unpaired
# rep1, 99929_2
Vienna_sle %>% filter(sample == "rep1" & readname == "99929_2") %>% 
  select(sample, readname,dG, piRNA,pilen, piRNAseq, targetlen,targetseq,piRNAmatch, targetmatch) %>% 
  write_csv(paste0(outdir,"CLASH_Ext_first_unpaired_3.csv"))
# rep1, 1482069_1
Vienna_sle %>% filter(sample == "rep1" & readname == "1482069_1") %>% 
  select(sample, readname,dG, piRNA,pilen, piRNAseq, targetlen,targetseq,piRNAmatch, targetmatch) %>% 
  write_csv(paste0(outdir,"CLASH_Ext_first_unpaired_4.csv"))

