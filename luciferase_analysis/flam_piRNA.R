.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

#######################################################################################################
# flam antisense piRNAs stats
#######################################################################################################
outdir <- "/media/hermione/piRNA_luc/flam/"

piRNA_seq <- read_tsv("/media/pericles/CLASH/piRNA/piRNA_CPM.tsv") %>%
  select(name,seq, pilen,CPM_avg)
sample_order <- c("flam0.6", "flam0.3","flam0.2","flam0.1","flam0.08")

# RNAplex
functional_Vienna2 <- read_csv(paste0("/media/hermione/piRNA_luc/figure2/","mdg1_construct_struct.csv"),
                               col_types = cols(dG = col_double())) %>% 
  mutate(sample = target %>% str_remove("mdg1_") %>% str_remove("anti_")) %>% 
  left_join(piRNA_seq %>% select(name,CPM_avg) %>% rename(piRNA = name,CPM = CPM_avg), by = "piRNA")

type_order <- c("perfect match","functional","not functional")
sample_order <- c("flam0.6", "flam0.3","flam0.2","flam0.1","flam0.08")
functional_Vienna2 %>% filter(str_detect(target,"flam")) %>% 
  mutate(target = target %>% str_replace("_","\\.")) %>% 
  select(pistruct, target,piRNA, CPM) %>% distinct() %>%
  group_by(pistruct, target) %>% summarise(sumCPM = sum(CPM)) %>% ungroup() %>%
  mutate(pistruct = pistruct %>% factor(levels = type_order), target = target %>% factor(levels = sample_order)) %>% 
  ggplot(aes(x = target, y = sumCPM, fill = pistruct)) +
  geom_col( color ="black", position = "stack") +
  labs(y = "total piRNA expression (CPM)")+
  scale_y_continuous(breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500)) +
  theme_minimal() + scale_fill_manual(values = c("maroon1","royalblue1","gray65"))+
  theme(axis.title.x = element_blank(), legend.title = element_blank(), panel.grid.minor.y = element_blank(),
        legend.position = c(0.7,0.8), axis.text.x = element_text(angle = 30, vjust = 1)) +
  ggsave(paste0(outdir,"piRNA_CPM_flam3.png"),  width =3, height =3, dpi = 300)+
  ggsave(paste0(outdir,"piRNA_CPM_flam3.pdf"),  width =3, height =3)


#######################################################################################################
# plasmid comparison
# piRNA is antisense, so 5' ssend -> sstart 3'.
#######################################################################################################
outdir <- "/media/hermione/piRNA_luc/plasmid/"

sample_order <- c("mdg1_anti", "piRNA134303", "piRNA268036")

plasmid_target <- read_tsv("/media/hermione/piRNA_luc/plasmid/piRNA_plasmid_target.tsv",
                                col_names = c("name","seq")) %>% 
  separate_rows(seq,sep ="") %>% filter(!(seq =="")) %>%
  group_by(name) %>% mutate(index = row_number() %>% as.integer(),
                            my_col = case_when(index <= 10 ~ "black",
                                               max(index) - index < 10  ~ "black",
                                               TRUE ~ "red")) %>% ungroup() %>% 
  mutate(label2 = paste0(seq, "\n",index))

# RNAplex
functional_Vienna2 <- read_csv(paste0("/media/hermione/piRNA_luc/figure2/","mdg1_construct_struct.csv"),
                               col_types = cols(dG = col_double())) %>% 
  mutate(sample = target %>% str_remove("mdg1_") %>% str_remove("anti_")) %>% 
  left_join(piRNA_seq %>% select(name,CPM_avg) %>% rename(piRNA = name,CPM = CPM_avg), by = "piRNA") %>% 
  filter(!(str_detect(target,"flam")))

type_order <- c("perfect match","functional","not functional")
sample_order <- c("mdg1_anti", "piRNA134303", "piRNA268036")
functional_Vienna2 %>% mutate(target = target %>% str_replace("pi_RNA","piRNA")) %>% 
  select(pistruct, target,piRNA, CPM) %>% distinct() %>%
  group_by(pistruct, target) %>% summarise(sumCPM = sum(CPM)) %>% ungroup() %>%
  mutate(pistruct = pistruct %>% factor(levels = type_order), target = target %>% factor(levels = sample_order)) %>% 
  ggplot(aes(x = target, y = sumCPM, fill = pistruct)) +
  geom_col( color ="black", position = "stack") +
  labs(y = "total piRNA expression (CPM)")+
  scale_y_continuous(breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500)) +
  theme_minimal() + scale_fill_manual(values = c("maroon1","royalblue1","gray65"))+
  theme(axis.title.x = element_blank(), legend.title = element_blank(), panel.grid.minor.y = element_blank(),
        legend.position = c(0.7,0.8), axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggsave(paste0(outdir,"piRNA_CPM_plasmid3.png"),  width =2.5, height =3, dpi = 300)+
  ggsave(paste0(outdir,"piRNA_CPM_plasmid3.pdf"),  width =2.5, height =3)

plasmid_visualise3 <- function(data, targetlen, samplename, height, outname, width = 5.5){
  plasmid_color <- plasmid_target %>% filter(name == samplename) 
  plasmid_region_min <- plasmid_color %>% filter(my_col == "red") %>% .$index %>% min()
  plasmid_region_max <- plasmid_color %>% filter(my_col == "red") %>% .$index %>% max()
  data %>% filter(target == samplename) %>% arrange(sstart,send) %>% 
    mutate(cumcount_en = cumsum(CPM), cumcount_st = cumcount_en- CPM) %>% 
    mutate(pistruct = pistruct %>% factor(levels = type_order)) %>%
    ggplot(aes(xmin = send-0.5,xmax = sstart-0.5,ymin = cumcount_st,ymax = cumcount_en, fill = pistruct))+
    annotate("rect", xmin = plasmid_region_min -0.5, 
             xmax = plasmid_region_max +0.5, ymin =-Inf, ymax = Inf, alpha = 0.3,fill = 'gold')+
    geom_rect() + theme_minimal()+
    scale_fill_manual(values = c("maroon1","royalblue1","gray65")) +
    scale_x_continuous(breaks=1:targetlen, labels = plasmid_color %>% .$seq, expand = c(0, 0), limits = c(0,targetlen+1))+
    labs(y = "CPM")+
    theme(axis.text.x = element_text(color = plasmid_color$my_col, family = "Courier", size = 6,face = "bold"),
          axis.title.x = element_blank(),legend.position = "bottom", panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 7), legend.title = element_blank())+
    ggsave(paste0(outdir,outname,"_cumsum3.png"), width =width, height =6, dpi = 300 )+
    ggsave(paste0(outdir,outname,"_cumsum3.pdf"), width =width, height =6 )
}
piRNA_plasmid3 <- functional_Vienna2 %>% mutate(target = target %>% str_replace("pi_RNA","piRNA")) %>% rowwise() %>% 
  mutate(test = list(str_locate(targetmatch,"\\)")) %>% .[[1]],
         test2 = list(str_locate(targetmatch,"\\)\\.+$")) %>% .[[1]] ) %>% ungroup() %>% 
  mutate(sstart = test[,1], send = test2[,1] +1,pilen = str_length(piRNAseq))

plasmid_visualise3(piRNA_plasmid3, targetlen = 47, 
                   samplename = "piRNA134303", height = 12,outname = "piRNA134303", width = 5)
plasmid_visualise3(piRNA_plasmid3, targetlen = 45, 
                   samplename = "piRNA268036", height = 12,outname = "piRNA268036", width = 5)
plasmid_visualise3(piRNA_plasmid3, targetlen = 48, 
                   samplename = "mdg1_anti", height = 12,outname = "mdg1_anti", width = 5)


functional_percentage2 <- function(dataname, outname, samplename) {
  dataname %>% filter(target == samplename) %>%  select(pistruct, CPM) %>%distinct() %>% group_by(pistruct) %>%
    summarise(n = sum(CPM)) %>% ungroup() %>%  mutate(prop = 100 *n/sum(n)) %>% 
    mutate(ypos = cumsum(prop)-0.5*prop)  %>% 
    mutate(pistruct = pistruct %>% factor(levels = type_order)) %>% 
    ggplot(aes(x="", y=prop, fill=pistruct)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0, direction = 1) +
    theme_void() + theme(legend.position = "bottom")+
    guides(fill= guide_legend(title = NULL))+ scale_fill_manual(values = c("maroon1","royalblue1","gray65")) + 
    geom_text(aes(x = 1.6, label = paste0(round(prop,1), "%")), position = position_stack(vjust = 0.5))+
    ggsave(paste0(outdir,"plasmid_",outname,"_pie2.png"), width =4, height =3, dpi = 300)+
    ggsave(paste0(outdir,"plasmid_",outname,"_pie2.pdf"), width =4, height =3)
}

functional_percentage2(dataname = piRNA_plasmid3 , outname = "mdg1_anti", samplename = "mdg1_anti") 

#######################################################################################################
# flam mapping of piRNA134303, piRNA268036 & mdg1_anti
# piRNA is sense.
#######################################################################################################

outdir <- "/media/hermione/piRNA_luc/flam_mapping/"

sample_order <- c("mdg1_anti", "piRNA134303", "piRNA268036")
origin_plasmid <- read_tsv("/media/hermione/piRNA_luc/flam_mapping/flam_mapping.txt",
                          col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  mutate(qseqid = qseqid %>% str_extract("piRNA\\d+")) %>%
  left_join(piRNA_seq %>% rename(qseqid = name), by = "qseqid") %>% 
  rename(CPM = CPM_avg) %>% filter(sstart < send) %>% filter(length == pilen) %>% 
  mutate(sample = sseqid %>% str_extract("\\w+_hap\\d+"), sseqid = sseqid %>% str_remove("\\w+_hap\\d+_"))

origin_plasmid_target <- read_tsv("/media/hermione/piRNA_luc/flam_mapping/flam_mapping_target.tsv",
                           col_names = c("name","seq")) %>% 
  mutate(sample = name %>% str_extract("\\w+_hap\\d+"), name = name %>% str_remove("\\w+_hap\\d+_")) %>% 
  separate_rows(seq,sep ="") %>% filter(!(seq =="")) %>%
  group_by(name,sample) %>% mutate(index = row_number() %>% as.integer(),
                            my_col = case_when(index < 32 & name != "mdg1" ~ "black",
                                               (max(index) - index < 30) & name != "mdg1" ~ "black",
                                               index < 52 & name == "mdg1" ~ "black",
                                               (max(index) - index < 50) & name == "mdg1" ~ "black",
                                               TRUE ~ "red")) %>% ungroup() %>% 
  mutate(label2 = paste0(seq, "\n",index))
flam_region_anno <- read_tsv(paste0("/media/hermione/piRNA_luc/process/", "region_info.tsv")) %>% 
  mutate(new_start = regstart - start+ 1.5, new_end = regend - start+ 1.5, new_label = paste0(sample,":",start, ":",end))
origin_visualise <- function(data, targetlen, samplename, height, outname, width = 5.5, hap){
  origin_region <- flam_region_anno %>% filter(name == samplename) %>% filter(sample == hap)
  origin_color <- origin_plasmid_target %>% filter(name == samplename) %>% filter(sample == hap)
  data %>% filter(sseqid == samplename) %>% filter(sample == hap) %>% arrange(sstart,send) %>% 
    mutate(cumcount_en = cumsum(CPM), cumcount_st = cumcount_en- CPM) %>% 
    ggplot(aes(xmin = sstart-0.5,xmax = send+0.5,ymin = cumcount_st,ymax = cumcount_en, fill = CPM))+
    annotate("rect", xmin = origin_region$new_start[1], 
             xmax = origin_region$new_end[1], ymin =-Inf, ymax = Inf, alpha = 0.3,fill = 'gray')+
    geom_rect() + theme_minimal()+ 
    scale_fill_gradient(low = "yellow",high = "blue",trans = "log1p",limits=c(0,1000), breaks = c(0,10,100,500,1000)) +
    scale_x_continuous(breaks=1:targetlen, labels = origin_color %>% .$seq, expand = c(0, 0), limits = c(0,targetlen+1))+
    labs(y = "CPM")+
    annotate("text", x = origin_region$new_end[1]+10, y = 25, label = paste0(origin_region$new_label[1]))+
    theme(axis.text.x = element_text(color = origin_color$my_col, family = "Courier", size = 6,face = "bold"),
          axis.title.x = element_blank(),legend.position = "bottom", panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 7),
          legend.key.width = unit(1.5, "cm"))+
    ggsave(paste0(outdir,hap, "_",outname,"_cumsum.png"), width =width, height =6, dpi = 300 )+
    ggsave(paste0(outdir,hap, "_",outname,"_cumsum.pdf"), width =width, height =6 )
}
origin_visualise(origin_plasmid, targetlen = 88,  hap = "flam_hap1",
                  samplename = "piRNA134303", height = 12,outname = "piRNA134303", width = 7)
origin_visualise(origin_plasmid, targetlen = 86, hap = "flam_hap1",
                  samplename = "piRNA268036", height = 12,outname = "piRNA268036", width = 7)
