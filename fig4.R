

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/CLASH/vienna/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(ggrepel)))
suppressMessages(suppressWarnings(require(ggseqlogo)))
suppressMessages(suppressWarnings(require(gt)))
ppi <- 600
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

outdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/fig4/"


#######################################################################################################
# insertion stat
# sup_figS4C
#######################################################################################################
TEgroup <- read_tsv("/media/pericles/TEfind/database/class_process.tsv") %>%
  unite(dis, name, sep = "-",col = "name") %>%
  select(name, family, subfamily)

pbsv_repmask2 <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_mapped_all.tsv")
pbsv_repmask3 <- pbsv_repmask2 %>%
  mutate(len = insend-insstart, score2 = score*len) %>%
  group_by(chr, start, end, name2, strand,TEname) %>%
  summarise(score3 = sum(score2), total_len = sum(len)) %>% ungroup() %>%
  mutate(new_score = score3/total_len) %>%
  group_by(chr, start,end) %>% filter(new_score == min(new_score)) %>% filter(total_len == max(total_len)) %>%  ungroup() %>%
  select(chr, start, end, TEname, new_score, strand,total_len) %>% distinct() %>%
  mutate(type = "pbsv")

rmblast_len <- read_tsv("/media/pericles/TEfind/rmblast/dm6_rmblast_mask.tsv") %>%
  rename(TEname = name) %>% mutate(TElen = TEend-TEstart)
rmblast_repmask <- read_tsv("/media/pericles/TEfind/vcfprocess/embl_pbsv_del_repeat.bed",
                            col_names = c("chr","start", "end", "TEname","new_score","strand")) %>%
  mutate(type = "rmblast") %>%
  left_join(rmblast_len %>% select(chr, start,end,strand, TEname,TElen), by = c("chr","start","end","TEname","strand"))

tree_repmask <- bind_rows(pbsv_repmask3 %>% rename(TElen = total_len), rmblast_repmask) %>%
  left_join(TEgroup %>% rename(TEname = name), by = "TEname") %>%
  mutate(TEname = paste0(TEname,"@",type))
tree_repmask %>% 
  write_tsv(paste0(outdir,"sup_fig4_inserted_TEs.tsv"))
level <- c( "pbsv","rmblast")
family_level <- c("DNA","RC","LTR","LINE","SINE","Unknown")

tree_repmask %>%filter(chr != "chrY") %>% 
  mutate(type = case_when(type == "rmblast" ~ "reference TE",
                          type == "pbsv" ~ "newly inserted TE",
                          TRUE ~ "non")) %>%
  mutate(family = family %>% factor(levels = family_level)) %>% 
  ggplot(aes(x = new_score, fill = family)) +
  geom_histogram(binwidth = 1)+
  theme_minimal()+ theme(legend.position = "bottom")+
  labs(x = "RepeatMasker Score", y = "Number of TEs") + guides(fill = guide_legend(title = "TE family"))+
  facet_wrap(~type, nrow = 1) +scale_fill_brewer(palette="Dark2", direction = 1) +
  ggsave(paste0(outdir, "sup_fig4C.png"), width =5, height =4, dpi = 300)+
  ggsave(paste0(outdir, "sup_fig4C.pdf"), width =5, height =4)

tree_repmask %>% filter(chr != "chrY") %>% 
  group_by(type, family) %>%
  summarise(n = n()) %>% ungroup() %>% 
  group_by(type) %>% 
  mutate(prop = 100 *n/sum(n)) %>% 
  mutate(ypos = cumsum(prop)-0.5*prop) %>% ungroup() %>% 
  mutate(family = family %>% factor(levels = family_level)) %>% 
  ggplot(aes(x="", y=prop, fill=family)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0, direction = 1) + facet_wrap(~type)+
  theme_void() +scale_fill_brewer(palette="Dark2", direction = 1) +
  guides(fill= "none")+
  geom_text(aes(label = if_else(prop > 3,paste0(round(prop,1), "%"),"")), position = position_stack(vjust = 0.5))+
  ggsave(paste0(outdir,"sup_fig4C_pie_family.png"), width =5, height =3, dpi = 300)+
  ggsave(paste0(outdir,"sup_fig4C_pie_family.pdf"), width =5, height =3)


#######################################################################################################
# comparison of other data
#######################################################################################################

cluster3 <- read_tsv("/media/pericles/TEfind/compare/cluster2_remove.bed",
                     col_names = c("chr","start","end","name","score","strand","group"))

cluster4 <- cluster3 %>% 
  mutate(sample= case_when(str_detect(name, "pbsv") ~ "pbsv",
                           TRUE ~ "tebreak")) %>% 
  select(sample, group) %>% distinct()

VennDiagram::venn.diagram(
  x = list(cluster4 %>% filter(sample =="pbsv") %>% pull(group),
           cluster4 %>% filter(sample =="tebreak") %>% pull(group)),
  category.names = c(paste0("pbsv (", cluster4 %>% filter(sample =="pbsv") %>% pull(sample) %>% length(),")"),
                     paste0("tebreak (", cluster4 %>% filter(sample =="tebreak") %>% pull(sample) %>% length(),")")),
  filename = paste0(outdir,'sup_TE_ins_venn.png'),
  output=TRUE,imagetype="png" , height = 480 , width = 480 , resolution = 300, 
  compression = "lzw", lwd = 1, col=c('#21908dff', '#ca0020'), fill = c( alpha('#21908dff',0.3), alpha('#ca0020',0.3)),
  cex = 0.5, fontfamily = "sans", cat.cex = 0.4, cat.default.pos = "outer",
  cat.pos = c(-27, 27), cat.dist = c(0.055, 0.055), cat.fontfamily = "sans", cat.col = c('#21908dff', '#ca0020'))
#######################################################################################################
# overview of insertion sites
#######################################################################################################
anno_cen_telo <- c("telomere","centromere","centromere","telomere","telomere","centromere","centromere","telomere","telomere","telomere", "telomere","centromere","telomere")
centromere_telomere <- read_tsv("/media/pericles/TEfind/database/centromere_telomere.txt",
                                col_names = c("chr","start","end","V4")) %>% 
  filter(chr %in% c("chr2L","chr2R","chr3L","chr3R","chr4", "chrX")) %>% 
  bind_cols(anno_cen_telo %>% enframe(value = "annotation")) %>% 
  dplyr::select(chr, start,end, annotation)

tree_repmask_out <- tree_repmask %>% separate(TEname, sep = "@",into =c("TEname","dis")) %>%
  group_by(type, family,TEname) %>% summarise(count = n()) %>% ungroup() %>%
  group_by(type) %>%
  arrange(desc(count)) %>%
  top_n(20) 
tree_repmask_out %>% write_csv(paste0(outdir,"sup_fig4B_TEcount.csv"))
tree_repmask_out %>%  gt() %>% gtsave(paste0(outdir,"sup_fig4B_TEcount.pdf"))
tree_repmask_out %>%  gt() %>% gtsave(paste0(outdir,"sup_fig4B_TEcount.png"))
tree_repmask_out %>%  gt() %>% gtsave(paste0(outdir,"sup_fig4B_TEcount.html"))

#######################################################################################################
# flam TE stats
#######################################################################################################
TEgroup <- read_tsv("/media/pericles/TEfind/database/class_process.tsv") %>%
  unite(dis, name, sep = "-",col = "name") %>%
  select(name, family, subfamily)
read_TE_file <- function(filename,dataname2) {
  datafile <-  read_tsv(paste0("/media/pericles/TEfind/flam/rmblast/",filename,"_TE.bed"), 
                            col_names = c("chr","start","end", "name", "quality", "strand")) %>% 
    separate(name, sep = "_", into = c("name", "TEstart", "TEend")) %>% 
    left_join(TEgroup, by = "name") %>% mutate(dataname = dataname2)
  return( datafile )
}
flam_hap1_TE <- read_TE_file(filename = "flam_hap1",dataname2 = "flam hap1")
flam_hap2_TE <- read_TE_file(filename = "flam_hap2",dataname2 =  "flam hap2")
l20A_hap1_TE <- read_TE_file(filename = "20A_hap1",dataname2 = "20A hap1")
l20A_hap2_TE <- read_TE_file(filename = "20A_hap2",dataname2 = "20A hap2")

family_order <- c("DNA", "LINE", "gypsy", "pao","copia","others")
datafile_order <- c("flam hap1", "flam hap2","20A hap1","20A hap2")
strand_order <- c("sense","antisense")

TE_all <- bind_rows(flam_hap1_TE,flam_hap2_TE,l20A_hap1_TE,l20A_hap2_TE)
TE_all %>% mutate(family2 = case_when(family == "DNA" | family == "RC" ~ "DNA",
                                        family == "LTR" & subfamily == "Copia" ~ "copia",
                                        family == "LTR" & subfamily == "Pao" ~ "pao",
                                        family == "LTR" & subfamily == "Gypsy" ~ "gypsy",
                                        family == "LTR" & subfamily == "Unknown" ~ "others",
                                        TRUE ~ family)) %>% 
  mutate(family2 = family2 %>% factor(levels = family_order), 
         dataname = dataname %>% factor(levels = datafile_order),
         strand = if_else(strand == "+", "sense","antisense") %>% factor(levels = strand_order),
          len = end-start +1) %>% 
  group_by(dataname,strand, family2,name) %>%  summarise(lencount = sum(len)/1000) %>% ungroup() %>% 
#  group_by(dataname) %>%  mutate(per = 100 * lencount / sum(lencount)) %>% ungroup() %>%
  filter(family2 != "others") %>% 
  ggplot(aes(y = fct_reorder(name, lencount) , x = lencount, fill = strand))+
  geom_col(position = "stack") + facet_grid(family2~dataname, scales = "free_y", space = "free_y")+
  theme_minimal()+ labs(x = "total length (kbp)")+
  theme(axis.title.y = element_blank(), legend.position = c(0.9,0.8), legend.title = element_blank()) +
  scale_fill_manual(values = c("red3","dodgerblue")) +
  ggsave(paste0(outdir,"sup_fig4_flam_TE.png"), width =7, height =8, dpi = 300)+
  ggsave(paste0(outdir,"sup_fig4_flam_TE.pdf"), width =7, height =8)

# percentage of cluster content

TE_all2_per <- TE_all %>% 
  mutate(dataname = dataname %>% factor(levels = datafile_order),
         strand = if_else(strand == "+", "sense","antisense") %>% factor(levels = strand_order),
         len = end-start +1) %>% 
  group_by(dataname,strand) %>%  summarise(lencount = sum(len)) %>% ungroup() %>% 
  group_by(dataname) %>%  mutate(prop = 100* lencount/sum(lencount)) %>% ungroup()

TE_all2_per %>% group_by(dataname) %>% mutate(ypos = cumsum(prop)-0.5*prop) %>% ungroup() %>% 
  ggplot(aes(x="", y=prop, fill=strand)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0, direction = 1) +
  theme_void() + theme(legend.position = "bottom")+
  guides(fill= guide_legend(title = NULL))+ scale_fill_manual(values = c("red3","dodgerblue")) + 
  geom_text(aes(x = 1.2, label = paste0(round(prop,1), "%")), position = position_stack(vjust = 0.5))+
  facet_wrap(~dataname, ncol = 1)+
  ggsave(paste0(outdir,"sup_fig4_flam_TE_per.png"), width =3, height =8, dpi = 300)+
  ggsave(paste0(outdir,"sup_fig4_flam_TE_per.pdf"), width =3, height =8)

# make genome browser
TE_all %>% mutate(strand = if_else(strand == "+", "sense","antisense"),
                  yposi_min = if_else(strand == "sense", 0,-1),
                  yposi_max = if_else(strand == "sense", 1,0),
                  start =start/1000,
                  end =end/1000) %>% 
  ggplot(aes(fill = strand))+
  geom_rect(aes(xmin = start, xmax = end, ymin = yposi_min,ymax =yposi_max))+
  theme_minimal() + 
  theme(axis.title.y = element_blank(), legend.position = c(0.9,0.8),
        axis.text.y = element_blank(),legend.title = element_blank(),
        panel.grid.major.y = element_blank()) +scale_y_continuous(breaks = c(-1,0,1))+
  labs(x = "distance from cluster's TSS (kbp)")+
  scale_fill_manual(values = c("dodgerblue","red3")) +
  facet_wrap(~chr, ncol = 1,strip.position ="left")+
  ggsave(paste0(outdir,"sup_fig4_flam_overall.png"), width =12, height =3.5, dpi = 300)
  
dependent_TE <- c("297","412","gypsy","gtwin","blood","17.6","mdg1","Tabor","I-element","HMS-Beagle2",
                  "Quasimodo", "Stalker","Stalker2","Stalker4","Idefix","mdg3","F-element","gypsy5","Juan")

TE_all2 <- TE_all %>% mutate(strand = if_else(strand == "+", "sense","antisense"),
                  yposi_min = if_else(strand == "sense", 0,-1),
                  yposi_max = if_else(strand == "sense", 1,0),
                  start =start/1000,
                  end =end/1000,
                  display = if_else((name %in% paste0("Dmel-",dependent_TE) & strand == "antisense") &
                                      (end-start) > 0.5, "display","no")) %>% 
 # filter(str_detect(chr,"flam")) %>% 
  filter(end < 100)
TE_all2 %>% ggplot(aes(fill = strand))+
  geom_rect(aes(xmin = start, xmax = end, ymin = yposi_min,ymax =yposi_max), color = "snow2")+
  geom_text_repel(data = subset(TE_all2, display == "display"),
                  aes(label = name,x = (start+end)/2,y=-0.5), size = 2.5,max.overlaps = Inf,
                  segment.color = "gray26", fontface = "bold", segment.size = 0.3,point.padding = 0.1,
                  box.padding = unit(0.5, "lines"), min.segment.length = 0,nudge_x = 0.1, nudge_y = 0.4)+
  theme_minimal() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.position = c(0.9,0.8), legend.title = element_blank(),
        panel.grid.major.y = element_blank()) +scale_y_continuous(breaks = c(-1,0,1))+
  labs(x = "distance from cluster's TSS (kbp)")+ 
  scale_fill_manual(values = c("deepskyblue","salmon1")) +
  facet_wrap(~chr, ncol = 1,strip.position ="left")+
  ggsave(paste0(outdir,"sup_fig4_flam_zoom.png"), width =12, height =4, dpi = 300)

#######################################################################################################
# flam piRNA stats
#######################################################################################################
genicpiRNA <- read_tsv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/piRNA/", "piRNA_CPM.tsv")) %>% 
  rename(piRNA = name,CPM = CPM_avg) %>% select(piRNA,CPM)


flam_piRNA_cal <- function(filename1 , bedsize = 100, binsize = 100) {
  piRNA_long_cov <- read_tsv(paste0("/media/pericles/TEfind/flam/new_piRNA/piRNA_all_",filename1,"_ctg.bed"),
                                  col_names = c("chr","start","end","name","qual","strand")) %>% 
    mutate(piRNA = name %>% str_extract("^piRNA\\d+")) %>% select(chr,start,end,piRNA,strand) %>% 
    left_join(genicpiRNA, by = "piRNA")
  
  piRNA_long_cov2 <- piRNA_long_cov %>% 
    mutate(start = binsize * (start %/% binsize), end = (start+binsize)) %>% 
    rename(bedstart = start, bedend = end) %>% 
    group_by(strand) %>% nest_legacy() %>%
    mutate(start = map(data,function(x) {x %>% group_by(bedstart) %>%
        summarise(plus=sum(CPM)) %>% ungroup() %>% rename(bed = bedstart)}),
        end = map(data,function(x) {x %>% group_by(bedend) %>%
            summarise(plus=sum(CPM)) %>% ungroup() %>% mutate(plus=-plus) %>% rename(bed = bedend)}),
        total = map2(start,end,function(x,y) {bind_rows(x,y) %>% group_by(bed) %>%
            summarise(count = sum(plus)) %>% ungroup() %>% arrange(bed) %>% mutate(cov = cumsum(count))})) %>%
    select(-data,-start,-end) %>% unnest_legacy()
  ggplot()+
    geom_rect(data= TE_all %>% mutate(strand = if_else(strand == "+", "sense","antisense"),
                                      yposi_min = if_else(strand == "sense", 0,-bedsize), yposi_max = if_else(strand == "sense", bedsize,0)) %>%
                filter(chr == filename1) %>% mutate(strand = if_else(strand == "sense","sense (TE)", "antisense (TE)")),
              aes(xmin = start, xmax = end, ymin = yposi_min,ymax =yposi_max, fill = strand), alpha = 0.4) + 
    geom_smooth(data = piRNA_long_cov2 %>% filter(strand == "+"),
                aes(x=bed,y=cov), color = "black", method = "loess",span = 0.1, alpha = 0.5)+
    geom_step(data = piRNA_long_cov2 %>% mutate(cov = if_else(strand == "-",-cov,cov)) %>% 
                mutate(strand = if_else(strand == "+", "sense mapped piRNA","antisense mapped piRNA")),
              aes(x=bed,y=cov,color = strand))+
    labs(x = "distance from cluster's TSS (bp)", y="piRNA expression (CPM)",fill="",color = "")+
    theme_minimal()+ scale_fill_manual(values = c("dodgerblue","red3")) + theme(legend.position = c(0.8,0.8)) +
    ggsave(paste0(outdir,"sup_fig4_long_piRNA_",filename1,".png"), width =8, height =4, dpi = 300)
  return(piRNA_long_cov2)
}

piRNA_long_flam_hap1 <- flam_piRNA_cal(filename1 = "flam_hap1",bedsize=1000, binsize = 100)
piRNA_long_flam_hap2 <- flam_piRNA_cal(filename1 = "flam_hap2",bedsize=1000, binsize = 100)
piRNA_long_20A_hap1 <- flam_piRNA_cal(filename1 = "20A_hap1",bedsize=100)
piRNA_long_20A_hap2 <- flam_piRNA_cal(filename1 = "20A_hap2",bedsize=100)


ggplot()+ geom_rect(data= TE_all %>% mutate(strand = if_else(strand == "+", "sense","antisense"),
                                    yposi_min = if_else(strand == "sense", 0,-500), yposi_max = if_else(strand == "sense", 500,0)) %>%
              filter(chr == "flam_hap1") %>% mutate(strand = if_else(strand == "sense","sense (TE)", "antisense (TE)")),
            aes(xmin = start, xmax = end, ymin = yposi_min,ymax =yposi_max, fill = strand), alpha = 0.4) + 
  geom_smooth(data = piRNA_long_flam_hap1 %>% filter(strand == "+"),
              aes(x=bed,y=cov), color = "black", method = "loess",span = 0.05, alpha = 0.5)+
  geom_step(data = piRNA_long_flam_hap1 %>% mutate(cov = if_else(strand == "-",-cov,cov)) %>% 
              mutate(strand = if_else(strand == "+", "sense mapped piRNA","antisense mapped piRNA")),
            aes(x=bed,y=cov,color = strand), alpha = 0.8)+
  labs(x = "distance from cluster's TSS (bp)", y="piRNA expression (CPM)",fill="",color = "")+ 
  coord_cartesian(ylim = c(-1000,2000),xlim = c(0,100000))+
  theme_minimal()+ scale_fill_manual(values = c("dodgerblue","red3")) + theme(legend.position = c(0.8,0.8)) +
  ggsave(paste0(outdir,"sup_fig4_long_piRNA_","flam_hap1","2.png"), width =8, height =4, dpi = 300)

ggplot()+ geom_rect(data= TE_all %>% mutate(strand = if_else(strand == "+", "sense","antisense"),
                                            yposi_min = if_else(strand == "sense", 0,-500), yposi_max = if_else(strand == "sense", 500,0)) %>%
                      filter(chr == "flam_hap2") %>% mutate(strand = if_else(strand == "sense","sense (TE)", "antisense (TE)")),
                    aes(xmin = start, xmax = end, ymin = yposi_min,ymax =yposi_max, fill = strand), alpha = 0.4) + 
  geom_smooth(data = piRNA_long_flam_hap2 %>% filter(strand == "+"),
              aes(x=bed,y=cov), color = "black", method = "loess",span = 0.05, alpha = 0.5)+
  geom_step(data = piRNA_long_flam_hap2 %>% mutate(cov = if_else(strand == "-",-cov,cov)) %>% 
              mutate(strand = if_else(strand == "+", "sense mapped piRNA","antisense mapped piRNA")),
            aes(x=bed,y=cov,color = strand), alpha = 0.8)+
  labs(x = "distance from cluster's TSS (bp)", y="piRNA expression (CPM)",fill="",color = "")+ 
  coord_cartesian(ylim = c(-1000,2000),xlim = c(0,100000))+
  theme_minimal()+ scale_fill_manual(values = c("dodgerblue","red3")) + theme(legend.position = c(0.8,0.8)) +
  ggsave(paste0(outdir,"sup_fig4_long_piRNA_","flam_hap2","2.png"), width =8, height =4, dpi = 300)


# y zoom

ggplot()+ geom_rect(data= TE_all %>% mutate(strand = if_else(strand == "+", "sense","antisense"),
                                            yposi_min = if_else(strand == "sense", 0,-500), yposi_max = if_else(strand == "sense", 500,0)) %>%
                      filter(chr == "flam_hap1") %>% mutate(strand = if_else(strand == "sense","sense (TE)", "antisense (TE)")),
                    aes(xmin = start, xmax = end, ymin = yposi_min,ymax =yposi_max, fill = strand), alpha = 0.4) + 
  geom_smooth(data = piRNA_long_flam_hap1 %>% filter(strand == "+"),
              aes(x=bed,y=cov), color = "black", method = "loess",span = 0.05, alpha = 0.5)+
  geom_step(data = piRNA_long_flam_hap1 %>% mutate(cov = if_else(strand == "-",-cov,cov)) %>% 
              mutate(strand = if_else(strand == "+", "sense mapped piRNA","antisense mapped piRNA")),
            aes(x=bed,y=cov,color = strand) , alpha = 0.8)+
  labs(x = "distance from cluster's TSS (bp)", y="piRNA expression (CPM)",fill="",color = "")+ 
  coord_cartesian(ylim = c(-1000,2000))+
  theme_minimal()+ scale_fill_manual(values = c("dodgerblue","red3")) + theme(legend.position = c(0.8,0.8)) +
  ggsave(paste0(outdir,"sup_fig4_long_piRNA_","flam_hap1","3.png"), width =10, height =4, dpi = 300)

ggplot()+ geom_rect(data= TE_all %>% mutate(strand = if_else(strand == "+", "sense","antisense"),
                                            yposi_min = if_else(strand == "sense", 0,-500), yposi_max = if_else(strand == "sense", 500,0)) %>%
                      filter(chr == "flam_hap2") %>% mutate(strand = if_else(strand == "sense","sense (TE)", "antisense (TE)")),
                    aes(xmin = start, xmax = end, ymin = yposi_min,ymax =yposi_max, fill = strand), alpha = 0.4) + 
  geom_smooth(data = piRNA_long_flam_hap2 %>% filter(strand == "+"),
              aes(x=bed,y=cov), color = "black", method = "loess",span = 0.05, alpha = 0.5)+
  geom_step(data = piRNA_long_flam_hap2 %>% mutate(cov = if_else(strand == "-",-cov,cov)) %>% 
              mutate(strand = if_else(strand == "+", "sense mapped piRNA","antisense mapped piRNA")),
            aes(x=bed,y=cov,color = strand), alpha = 0.8)+
  labs(x = "distance from cluster's TSS (bp)", y="piRNA expression (CPM)",fill="",color = "")+ 
  coord_cartesian(ylim = c(-1000,2000))+
  theme_minimal()+ scale_fill_manual(values = c("dodgerblue","red3")) + theme(legend.position = c(0.8,0.8)) +
  ggsave(paste0(outdir,"sup_fig4_long_piRNA_","flam_hap2","3.png"), width =10, height =4, dpi = 300)




#######################################################################################################
# piRNA expression profiles
#######################################################################################################

tx2 <- read_tsv(paste0("/media/pericles/CLASH/piRNA/","piRNA_CPM.tsv"))
# >piRNA134303
# TCCCTACAAGCTGTATGACCAAACCAT
# >piRNA268036
# TTTGGTCACTCTGACGCTCGCGTTG
# >mdg1_repeat
# TTTTTTATTTGTGGTTTTTTATTTGTGG

mdg1_functional <- read_tsv("/media/hermione/piRNA_luc/flam_mapping/flam_mapping.txt",
                            col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                          "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  mutate(qseqid = qseqid %>% str_extract("piRNA\\d+")) %>%
  left_join(tx2 %>% select(name, pilen) %>% rename(qseqid = name), by = "qseqid") %>% 
  mutate(sample = sseqid %>% str_extract("\\w+_hap\\d+"), sseqid = sseqid %>% str_remove("\\w+_hap\\d+_")) %>% 
  mutate(type = case_when((pilen == length & gapopen == 0) & mismatch==0 ~ "functional",
                          (qstart %in% c(1,2) & qend >= 20) & (gapopen == 0 & mismatch==0) ~ "functional",
                          TRUE ~ "not functional")) %>% 
  filter(sseqid == "mdg1") %>% 
  filter(type == "functional")

piRNAfromorder <- c("rmst","flam","20A","TE","genic","others")
piRNAgroup <- read_tsv(paste0("/media/pericles/CLASH/piRNA/", "anno/piRNA_all_annotation.tsv"))
tx4 <- tx2 %>%
  left_join(piRNAgroup, by = "name") %>% 
  mutate(type = if_else(type =="piRNA", "others",type)) %>% 
  mutate(type = type %>% str_remove("piRNA") %>% str_replace("TE2","TE")) 


cumsum_piRNA <- tx4 %>% mutate(type_raw = type) %>% 
  left_join(mdg1_functional %>% rename(name = qseqid) %>%
                                    select(name, sseqid) %>% distinct(), by = "name") %>% 
  mutate(type2 = case_when(seq == "TCCCTACAAGCTGTATGACCAAACCAT" ~ "piRNA134303",
                           seq == "TTTGGTCACTCTGACGCTCGCGTTG" ~ "piRNA268036",
                           sseqid == "mdg1" ~ "mdg1repeat",
                           type_raw == "flam" ~ "flam piRNA",
                           TRUE ~ "other piRNA"),
         type = case_when(type2 == "mdg1repeat" ~ "mdg1 repeat", 
                          type2 == "piRNA134303" ~ "flam piRNA", 
                          type2 == "piRNA268036" ~ "flam piRNA", 
                          type2 == "flam piRNA" ~ "flam piRNA", 
                          TRUE ~ "other piRNA")) %>% 
  group_by(type) %>% arrange(CPM_avg) %>% mutate(order = row_number() / n()) %>%
  ungroup()%>% mutate(logCPM = log10(CPM_avg))%>% group_by(type) %>% 
  mutate(allcount = n()) %>% ungroup() %>% 
  mutate(type3 = paste0(type," (n=",allcount,")"))


cumsum_piRNA  %>% 
  ggplot(aes(x = logCPM, y = order, color = type3))+
  annotate('rect',xmin = log10(5),xmax = Inf,ymin = 0,ymax = 1,alpha = 0.1, fill = "orangered")+
  annotate('rect',xmin = 0,xmax =log10(5),ymin = 0,ymax = 1,alpha = 0.1, fill = "aquamarine")+
  geom_step() + 
  ggrepel::geom_text_repel(data = subset(cumsum_piRNA %>% filter(type2 %in% c("piRNA134303", "piRNA268036"))),
                           aes(label =type2), size = 3, color = "black",  point.padding = 0.1,
                           box.padding = unit(4, "lines"), min.segment.length = 0)+
 # annotate("text",x = 2,y = 0.5,label = paste0("mdg1 repeat = ", cumsum_piRNA %>% filter(type2 == "mdg1repeat") %>% .$name %>% length()))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_minimal() +theme(legend.position = c(0.5,0.2), legend.title = element_blank())+
  labs(x = "log10 piRNA expression (CPM)", y = "Cumulative fraction") +
  ggsave(paste0(outdir,"sup_fig6_piRNA_express_cumsum.png"), width =3.5, height =3.5, dpi = 300)+
  ggsave(paste0(outdir,"sup_fig6_piRNA_express_cumsum.pdf"), width =3.5, height =3.5)
  
#######################################################################################################
# piRNA coming stats
#######################################################################################################
tx4 <- read_tsv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/piRNA/", "piRNA_CPM.tsv"))


# CPM
tx4 %>% 
  select(name,CPM_avg, type3) %>% distinct() %>%
  group_by(type3) %>% summarise(CPM = sum(CPM_avg)) %>% ungroup() %>%
  mutate(prop = CPM * 100/sum(CPM)) %>%
  mutate(ypos = cumsum(prop)-0.5*prop) %>%
  mutate(type3 = type3 %>% factor(levels = piRNAfromorder)) %>% 
  ggplot(aes(x="", y=prop, fill=type3)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0, direction = 1) +
  theme_void() +
  guides(fill= guide_legend(title = NULL))+
  geom_text(aes(label = if_else(prop > 3,paste0(round(prop,1), "%"),"")), position = position_stack(vjust = 0.5), size = 2.5)+
  scale_fill_brewer(palette="Set3", direction = -1) +
  ggsave(paste0(outdir,"sup_fig4_piRNA_from_pie_CPM.png"), width =4, height =4, dpi = 300)+
  ggsave(paste0(outdir,"sup_fig4_piRNA_from_pie_CPM.pdf"), width =4, height =4)


# main figures
#######################################################################################################
# Fig4B
#######################################################################################################
# CAGE
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
         logFC = log2(siPIWI/siEGFP)) %>% select(-contains("rep"))




# clustering data analysis
# CAGE

cluster_CLASH2 <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/fig2/", "clash_sle_cluster3.csv"), col_types = cols(dG = col_double())) %>% 
  filter(targettype%in% c("mRNA","ensemblTE")) %>% select(sample,readname,piRNA, cluster,targettype) %>% drop_na() %>% 
  left_join(Vienna_select, by = c("sample","readname","piRNA","targettype")) %>% drop_na() %>% distinct() %>% 
  left_join(CAGE_table %>% rename(target = genename), by = "target") %>% drop_na() %>% 
  group_by(target) %>% mutate(totalcount = n()) %>% ungroup() %>% filter(totalcount >= 2) %>% 
  group_by(sample,readname) %>% mutate(remover = n()) %>% ungroup() %>% filter(remover < 3) %>% 
  mutate(matetype = case_when(!(str_detect(target,"FBgn")) ~ "TE",TRUE ~ "mRNA")) %>% 
  group_by(target,cluster,matetype) %>% count() %>% ungroup() %>% 
  group_by(target) %>% slice(which.max(n)) %>% ungroup() %>% distinct() %>% 
  mutate(matetype2 = paste0(matetype, " (",cluster,")")) %>% select(-n)

mRNA_table3 <- CAGE_table %>% left_join(cluster_CLASH2 %>% rename(genename = target), by = "genename") %>% 
  mutate(matetype2 = case_when(!(is.na(matetype2)) ~ matetype2,
                              !(str_detect(genename,"FBgn")) ~ "TE (not seen)",
                              TRUE ~ "mRNA (not seen)"),
         type = if_else(str_detect(matetype2, "TE"), "TE", "mRNA")) 

mRNA_table4 <- mRNA_table3 %>% filter(!(is.na(logFC) | is.infinite(logFC))) %>%
  group_by(type,matetype2) %>% arrange(matetype2,logFC) %>% mutate(order = row_number() / n()) %>% ungroup() %>%
  mutate(change = if_else(logFC >= 1,"plus_change",if_else(logFC <= -1,"minus_change","stay")))

mRNAtypeorder <- c("TE (cluster1)","TE (cluster2)","TE (cluster3)","TE (not seen)",
                   "mRNA (cluster1)","mRNA (cluster2)","mRNA (cluster3)","mRNA (not seen)")
mRNA_table4 %>% mutate(matetype2 = matetype2 %>% factor(levels = mRNAtypeorder)) %>% 
  ggplot(aes(x = logFC, y = order, color = matetype2))+
  geom_step(alpha = 0.6,size =1) + scale_color_brewer(palette = "Set1")+
  theme_minimal()+ labs(x = "log2FC (CAGE)", y = "Cumulative Fraction")+
  theme(legend.position = c(0.2,0.6)) + guides(color=guide_legend(title="mRNA type \n(CLASH cluster)", size = 1.5))+
  ggsave(paste0(outdir,"scatter_cumsum","_cluster.png"), width =4, height =4, dpi = 300)+
  ggsave(paste0(outdir,"scatter_cumsum","_cluster.pdf"), width =4, height =4)


#######################################################################################################
# changed CAGE site, chimera tables
#######################################################################################################



CAGE_table2 <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "CAGE_edgeR.csv"))

close_bed <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_closest.bed",
                      col_names = c("chr1","start1","end1","ins_name","score1","strand1","chr2","start2","end2","gene","score2","strand2")) %>%
  mutate(state = "close")
inside_bed <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_inside.bed",
                       col_names = c("chr1","start1","end1","ins_name","score1","strand1","chr2","start2","end2","gene","score2","strand2")) %>%
  mutate(state = "inside")
inside_gene <- read_tsv("/media/pericles/TEfind/newTE/OSC_inside.bed",
                       col_names = c("chr1","start1","end1","ins_name","score1","strand1","chr2","start2","end2","gene","score2","strand2")) %>% 
  mutate(length = end1-start1) %>% filter(length != 0)

pbsv_repmask <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/pbsv_vcf_ins_read.fasta.out.gff",
                         col_names = c("name2","source", "dis1","insstart","insend", "score","strand","dis2","name"),skip = 3) %>%
  separate(name,into = c("dis3","TEname","TEstart","TEend"), sep = " ") %>%
  mutate(TEname = TEname %>% str_remove_all('"|Motif:'))
# direction of gene
gene_direction <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_tss.tsv") %>% 
  select(genesymbol,direction)


#change_CAGE <- CAGE_table2 %>% 
#  filter(str_detect(name,"FBgn")) %>% filter(fdr2 == "change" & logFC > 0) %>% 
#  left_join(close_bed %>% select(gene,ins_name) %>% rename(closeins = ins_name,name = gene), by = "name") %>% 
#  left_join(inside_bed %>% select(gene,ins_name) %>% rename(insideins = ins_name,name = gene), by = "name") %>% 
#  left_join(inside_gene %>% select(gene,ins_name) %>% rename(rmblastins = ins_name,name = gene) %>% distinct(), by = "name") %>% 
#  left_join(pbsv_repmask %>% select(name2, TEname) %>% rename(closeins = name2,closeTE = TEname)%>% distinct(), by = "closeins") %>% 
#  left_join(pbsv_repmask %>% select(name2, TEname) %>% rename(insideins = name2,insideTE = TEname)%>% distinct(), by = "insideins")


close_CAGE <- CAGE_table2 %>% filter(str_detect(name,"FBgn")) %>% filter(fdr2 == "change" & logFC > 0) %>% 
  left_join(close_bed %>% select(gene,ins_name) %>% rename(ins = ins_name,name = gene), by = "name") %>% 
  left_join(pbsv_repmask %>% select(name2, TEname,strand) %>% rename(ins = name2,TE = TEname,strand1 = strand)%>% distinct(), by = "ins") %>% 
  select(name, genesymbol,logFC,logCPM,ins,TE,FDR,strand1) %>% mutate(type = "close") %>% drop_na()
inside_CAGE <- CAGE_table2 %>% filter(str_detect(name,"FBgn")) %>% filter(fdr2 == "change" & logFC > 0) %>% 
  left_join(inside_bed %>% select(gene,ins_name) %>% rename(ins = ins_name,name = gene), by = "name") %>% 
  left_join(pbsv_repmask %>% select(name2, TEname,strand) %>% rename(ins = name2,TE = TEname,strand1 = strand)%>% distinct(), by = "ins") %>% 
  select(name, genesymbol,logFC,logCPM,ins,TE,FDR,strand1) %>% mutate(type = "inside") %>% drop_na()
rmblast_CAGE <- CAGE_table2 %>% filter(str_detect(name,"FBgn")) %>% filter(fdr2 == "change" & logFC > 0) %>% 
  left_join(inside_gene %>% select(gene,ins_name,strand1) %>% rename(ins = ins_name,name = gene) %>% distinct(), by = "name") %>% 
  select(name, genesymbol,logFC,logCPM,ins,FDR,strand1) %>% mutate(TE = ins,type = "rmblast")%>% drop_na()

change_CAGE <- bind_rows(close_CAGE, inside_CAGE, rmblast_CAGE)


Isoseq_chimera <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/","Isoseq_chimera_list.csv")) 
CAGE_chimera <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/","CAGE_all_chimera_list.csv")) %>% 
  pivot_wider(names_from = "sample",values_from = "count", values_fill =list(count=0)) %>% 
  mutate(siEGFP = siEGFPrep1 + siEGFPrep2+siEGFPrep3,
         siPIWI = siPIWIrep1 + siPIWIrep2+siPIWIrep3) %>% 
  select(-contains("rep")) %>% 
  filter(symbol != "notannotated") %>% filter(rnatype != "rmstRNA") %>% 
  group_by(symbol,  gene_name, TEtype) %>% summarise(siEGFP = sum(siEGFP),siPIWI = sum(siPIWI)) %>% ungroup() %>% 
  mutate(FC = log2((siPIWI+1)/(siEGFP +1)))

# CAGE chimera
change_CAGE2 <- change_CAGE %>% 
  left_join(CAGE_chimera %>% rename(name = gene_name), by = "name") %>% 
  filter(TE == TEtype)

All_fusion <- bind_rows(
  CAGE_table2 %>% filter(str_detect(name,"FBgn")) %>% filter(fdr2 == "change" & logFC > 0) %>% 
    select(name, genesymbol,logFC,logCPM,FDR) %>% anti_join(change_CAGE, by = "name"),
  change_CAGE %>% anti_join(change_CAGE2, by = "name"),
  change_CAGE2
) %>% left_join(gene_direction, by = "genesymbol")
All_fusion %>% write_csv(paste0(outdir,"Insertion_chimera_table.csv"))
 
# Isoseq chimera
Isoseq_chimera2 <- Isoseq_chimera %>% filter(state == "TE-genome" | state == "genome-TE") %>% 
  filter(!(str_detect(first,"novel") | str_detect(second,"novel"))) %>% 
  mutate(gene_name = if_else(str_detect(first, "FBgn"), first, second),
         symbol = if_else(str_detect(first, "FBgn"), firstsymbol, secondsymbol),
         TEtype = if_else(str_detect(first, "FBgn"), second, first)) %>% 
  select(gene_name, symbol,TEtype, KD,count_fl) %>% 
  group_by(gene_name, symbol,TEtype, KD) %>% summarise(count = sum(count_fl)) %>% ungroup() %>% 
  pivot_wider(names_from = "KD",values_from = "count")

change_Isoseq2 <- change_CAGE %>% 
  left_join(Isoseq_chimera2 %>% rename(name = gene_name), by = "name") %>% 
  filter(TE == TEtype)

All_fusion2 <- bind_rows(
  CAGE_table2 %>% filter(str_detect(name,"FBgn")) %>% filter(fdr2 == "change" & logFC > 0) %>% 
    select(name, genesymbol,logFC,logCPM,FDR) %>% anti_join(change_CAGE, by = "name"),
  change_CAGE %>% anti_join(change_Isoseq2, by = "name"),
  change_Isoseq2
) %>% left_join(gene_direction, by = "genesymbol") %>% distinct()
All_fusion2 %>% write_csv(paste0(outdir,"Insertion_chimera_table_Isoseq.csv"))
 
#######################################################################################################
# Previous data comparison
#######################################################################################################

anotation1 <- read_tsv("/media/pericles/CLASH/database/gtf/dm6_use.tsv") %>%
  select(type, name, genesymbol) %>%
  filter(type %in% c("rRNA", "tRNA", "snRNA", "snoRNA", "mRNA", "pseudogene", "ncRNA", "pre_miRNA", "miRNA")) %>%
  distinct() %>% select(name, genesymbol)


CAGE_table_raw <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "CAGE_Matrix.csv")) %>% 
  filter(str_detect(genename,"FBgn")) %>% 
  select(genename,genesymbol, sample,count) %>% 
  group_by(genename,genesymbol, sample) %>% summarise(count = sum(count)) %>% ungroup() %>% arrange(genename,sample) %>% 
  pivot_wider(names_from = "sample",values_from = "count" , values_fill =list(count=0)) %>% 
  mutate(siEGFP = siEGFP_rep2 +siEGFP_rep3  + siEGFP_rep1,
         siPIWI = siPIWI_rep1 + siPIWI_rep2 +siPIWI_rep3) %>% select(-contains("rep"))

Isoseq_table <- read_csv(paste0("/home/qqprosperodd/Desktop/Piwi_article_analysis/RNA/", "Isoseq_Matrix.csv"))


close_CAGE2 <- close_bed %>% select(gene,ins_name) %>% rename(ins = ins_name,name = gene) %>% 
  left_join(pbsv_repmask %>% select(name2, TEname,strand) %>% rename(ins = name2,TE = TEname,strand1 = strand)%>% distinct(), by = "ins") %>% 
  select(name, ins,TE,strand1) %>% mutate(type = "close") %>% drop_na()
inside_CAGE2 <- inside_bed %>% select(gene,ins_name) %>% rename(ins = ins_name,name = gene) %>% 
  left_join(pbsv_repmask %>% select(name2, TEname,strand) %>% rename(ins = name2,TE = TEname,strand1 = strand)%>% distinct(), by = "ins") %>% 
  select(name, ins,TE,strand1) %>% mutate(type = "inside") %>% drop_na()
rmblast_CAGE2 <- inside_gene %>% select(gene,ins_name,strand1) %>% rename(ins = ins_name,name = gene) %>% distinct() %>% 
  select(name, ins,strand1) %>% mutate(TE = ins,type = "rmblast")%>% drop_na()

change_CAGE3 <- bind_rows(close_CAGE2, inside_CAGE2, rmblast_CAGE2)




change_CAGE4 <- change_CAGE3 %>% 
  left_join(CAGE_chimera %>% rename(name = gene_name,CAGE_chimera_EGFP = siEGFP, CAGE_chimera_PIWI = siPIWI), by = "name") #%>% filter(TE == TEtype)

change_Isoseq4 <- change_CAGE3 %>% 
  left_join(Isoseq_chimera2 %>% mutate(Isoseq_TEtype = TEtype) %>% 
    rename(name = gene_name, Isoseq_symbol = symbol, Isoseq_chimera_EGFP = EGFP, Isoseq_chimera_PIWI = PIWI), by = "name") #%>% filter(TE == TEtype)

All_fusion_map <-  full_join(change_CAGE4,change_Isoseq4,
                            by = c("name",  "ins", "TE",  "strand1", "type", "TEtype")) %>% distinct()
change_CAGE5 <- change_CAGE3 %>% anti_join(All_fusion_map, by = "name")

All_fusion_3 <- anotation1  %>% 
  left_join(CAGE_table_raw %>% rename(name = genename), by = c("name","genesymbol")) %>% 
  left_join(bind_rows(change_CAGE5, All_fusion_map ), by = c("name")) %>% 
  left_join(Isoseq_table %>% select(associated_gene ,fl_EGFP, fl_PIWI) %>% rename(name = associated_gene), by = "name") %>% distinct() %>% 
  left_join(gene_direction, by = "genesymbol") %>% distinct()


strandorder <- c("sense", "antisense")


test <- All_fusion_3 %>% filter(TE == TEtype)


All_fusion_3 %>% filter(TE == TEtype) %>% distinct() %>% 
  mutate(strandness = if_else(strand1 == direction, "sense","antisense") %>% factor(levels = c("sense","antisense"))) %>% 
  group_by(genesymbol,TE, siEGFP,siPIWI) %>% 
  slice_max(n=1, desc(strandness), with_ties = FALSE) %>% ungroup() %>% 
  group_by(strandness) %>% count() %>% ungroup() %>% 
  mutate(total = sum(n), prop = 100 * n/total) %>% 
  mutate(ypos = cumsum(prop)-0.5*prop) %>%
  mutate(strandness = strandness %>% factor(levels = strandorder)) %>% 
  ggplot(aes(x="", y=prop, fill=strandness)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0, direction = 1) +
  theme_void() +   scale_fill_manual(values = c("deepskyblue","salmon1")) +
  guides(fill= guide_legend(title = NULL))+
  geom_text(aes(label = if_else(prop > 3,paste0(round(prop,1), "%"),"")), position = position_stack(vjust = 0.5), size = 2.5)+
  ggsave(paste0(outdir,"insert_strandness_CAGE.png"), width =4, height =4, dpi = 300)+
  ggsave(paste0(outdir,"insert_strandness_CAGE.pdf"), width =4, height =4)

All_fusion_3 %>% filter(TE == Isoseq_TEtype) %>% distinct() %>% 
  mutate(strandness = if_else(strand1 == direction, "sense","antisense") %>% factor(levels = c("sense","antisense"))) %>% 
  group_by(genesymbol,TE, fl_EGFP,fl_PIWI) %>% 
  slice_max(n=1, desc(strandness), with_ties = FALSE) %>% ungroup() %>% 
  group_by(strandness) %>% count() %>% ungroup() %>% 
  mutate(total = sum(n), prop = 100 * n/total) %>% 
  mutate(ypos = cumsum(prop)-0.5*prop) %>%
  mutate(strandness = strandness %>% factor(levels = strandorder)) %>% 
  ggplot(aes(x="", y=prop, fill=strandness)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0, direction = 1) +
  theme_void() +   scale_fill_manual(values = c("deepskyblue","salmon1")) +
  guides(fill= guide_legend(title = NULL))+
  geom_text(aes(label = if_else(prop > 3,paste0(round(prop,1), "%"),"")), position = position_stack(vjust = 0.5), size = 2.5)+
  ggsave(paste0(outdir,"insert_strandness_Isoseq.png"), width =4, height =4, dpi = 300)+
  ggsave(paste0(outdir,"insert_strandness_Isoseq.pdf"), width =4, height =4)

# previous DEG read

previous_DEG <- readxl::read_xlsx("/home/qqprosperodd/Desktop/Piwi_article_analysis/DEG_previous.xlsx")

All_fusion3_table <- previous_DEG %>% select(gene_name) %>% distinct() %>% rename(genesymbol = gene_name) %>% 
  left_join(All_fusion_3, by = "genesymbol")


All_fusion3_table %>% 
  write_csv(paste0(outdir,"DEG_table_raw.csv"))









