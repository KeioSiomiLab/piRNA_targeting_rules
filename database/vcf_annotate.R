#######################################################################################################
# annotate TE
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gt)))
suppressMessages(suppressWarnings(require(ggpointdensity)))
suppressMessages(suppressWarnings(require(viridis)))
ppi <- 300
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
TE <- read_tsv("/media/pericles/TEfind/newTE/refTE.tsv") %>% 
  mutate(length = str_length(seq),
         start = 0) %>% rename(TEname = key) %>%
  filter(str_detect(TEname, "Dmel-"))%>% mutate(TEname = TEname %>% str_remove("Dmel-"))
ref_anno <- read_tsv("/media/pericles/TEfind/database/annotateTE/TE_bed6_full.bed",
                     col_names = c("TEname","querystart","queryend","annotation","score","strand"))
ref_anno2 <- ref_anno %>% 
  group_by(TEname) %>% 
  summarise(refannotation = paste0(annotation,collapse = "&"),
            total = n()) %>% ungroup()


ref_anno %>% filter(!(annotation %in% c("start_codon","polyA_signal_sequence","TATA_box", "intron",
                                        "five_prime_UTR","three_prime_UTR","primer_binding_site", "transcription_start_site",
                                        "polyA_site","polyA_sequence"))) %>% 
  filter(str_detect(TEname, "Dmel-")) %>% mutate(TEname = TEname %>% str_remove("Dmel-")) %>% 
  mutate(annotation = annotation %>% str_remove("CDS_") %>% str_replace("terminal_inverted_repeat","TIR") %>% 
           str_replace("three_prime_LTR","3'LTR") %>% str_replace("five_prime_LTR","5'LTR") %>% 
           str_replace("direct_repeat","DIR") %>% str_replace("long_terminal_repeat","LTR")) %>% 
  ggplot(aes(y = TEname))+
  geom_segment(data = TE ,
               aes(x=start, xend = length,yend = TEname), size = 4,color = "gray", alpha = 0.4)+
  geom_segment(aes(x = querystart,xend = queryend, yend = TEname, color = annotation), size = 9, alpha = 0.7)+
  ggrepel::geom_text_repel(aes(x = (querystart + queryend)/2, label = paste0(annotation, "\n",querystart,"-",queryend)), 
                           color = "black",size = 3.2, lineheight = 1, nudge_x = 0, nudge_y = 0,max.overlaps = Inf, 
                           direction = "x", point.size = NA, force_pull = 0,
                           box.padding = 0)+
  theme_minimal() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000))+
  guides(color = FALSE)+
  ggsave(paste0("/media/pericles/TEfind/database/TE_annotation.png"), width =15, height =40, dpi = 300)




TEgroup <- read_tsv("/media/pericles/TEfind/database/class_process.tsv") %>% 
  unite(dis, name, sep = "-",col = "name") %>% 
  select(name, family, subfamily) %>% rename(TEname =name)
pbsv_full <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_locus.bed",
                      col_names = c("TEname","start","end","name"))

pbsv_annotate <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_annotate.bed",
                          col_names = c("queryTE","querystart","queryend","annotation","dis","querystrand",
                                        "TEname","start","end","name","overlap")) %>% 
  mutate(type = "pbsv")
embl_full <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_rmblast_locus.bed",
                      col_names = c("TEname","start","end","name"))

embl_annotate <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_rmblast_annotate.bed",
                          col_names = c("queryTE","querystart","queryend","annotation","dis","querystrand",
                                        "TEname","start","end","name","overlap")) %>% 
  mutate(type = "embl")

annotate <- bind_rows(pbsv_annotate,embl_annotate)

annotate2 <- annotate %>% 
  mutate(rate = overlap /(queryend-querystart)) %>% 
  filter(rate > 0.2 ) %>% 
  mutate(state = case_when(rate == 1 ~ "full",
                           rate > 0.9 ~ "almost_full",
                           rate > 0.5 ~ "fragile",
                           TRUE ~ "dust"))


annotate3 <- bind_rows(pbsv_full,embl_full) %>% left_join(annotate2,by = c("TEname","start","end","name")) %>% 
  group_by(name,TEname, start,end) %>% 
  summarise(annotation = paste0(annotation,collapse = "&"),
            state = paste0(state,collapse = "&"),
            count = n()) %>% 
  ungroup() %>% 
  mutate(LTR = if_else(str_detect(annotation,"LTR"),"LTR","no"),
         count = if_else(annotation=="NA",0L,count)) %>% 
  left_join(ref_anno2,by = "TEname") %>% 
  mutate(totalness = count/total) %>% 
  left_join(TEgroup, by = "TEname")
  

annotate3 %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/OSC_repeat_annotation.tsv")
#######################################################################################################
# make allTE file
#######################################################################################################
# must edit!
suppressMessages(suppressWarnings(require(tidyverse)))
TElist <- read_tsv("/media/pericles/TEfind/database/TE_dm.tsv",
                   col_names = c("name","seq")) %>% 
  mutate(key = name) %>%   select(key,name,seq)
genome_TE <- read_tsv("/media/pericles/TEfind/rmblast/TE_genome.tsv",
                      col_names = c("name2","seq")) %>% 
  separate(name2, sep = "\\(", into = c("name2","strand")) %>% 
  mutate(name2 = paste0(name2 %>% str_replace_all("-|\\:", "\\|"), "|",strand %>% str_replace_all("\\)", "\\|"), "rmblast")) %>% 
  select(name2, seq)

del_repeat <- read_tsv("/media/pericles/TEfind/vcfprocess/embl_pbsv_del_repeat.bed",
                       col_names = c("chr","start","end","name","score","strand"))

rmblast2 <- read_tsv("/media/pericles/TEfind/rmblast/dm6_rmblast_mask.tsv") %>% 
  mutate(name2 = paste0(chr, "|",start,"|",end,"|",strand,"|rmblast"))
rmblast3 <- rmblast2 %>% select(name2,name) %>% left_join(genome_TE,by = "name2") %>% 
  rename(key = name,name = name2) %>% 
  select(key, name, seq) %>% drop_na()

pbsv_ins <- read_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins_read.tsv",
                     col_names = c("name2","seq"))
pbsv_out  <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_anno.tsv") 

pbsv_TE <- pbsv_out %>% 
  left_join(pbsv_ins, by= "name2") %>% 
  mutate(seq = if_else(strand =="-",chartr("ATGCatgc","TACGtacg",seq) %>% stringi::stri_reverse(),seq)) %>% 
  mutate(name2 = paste0(chr, "|",start,"|",end,"|",strand,"|pbsv"),len = str_length(seq)) %>%
  group_by(name2) %>% arrange(desc(len)) %>% mutate(fil = row_number()) %>% ungroup() %>% 
  filter(fil ==1) %>% rename(key = TEname, name = name2) %>% select(key, name,seq)

tes <- bind_rows(TElist,pbsv_TE,rmblast3) %>% mutate(seq = seq %>% str_to_upper()) %>% 
  mutate(index = row_number(),name = paste0(">", name)) %>% 
  gather(seq,name,key = "type",value = "fasta") %>% 
  arrange(index, type) %>% select(key,fasta)
tes %>% select(fasta) %>% write_tsv("/media/pericles/TEfind/newTE/allTE.fasta",col_names = FALSE)



tes4 <- bind_rows(TElist,pbsv_TE,rmblast3) %>% mutate(seq = seq %>% str_to_upper(), len = str_length(seq)) %>% 
  group_by(key) %>% mutate(reflen = str_length(seq) %>% .[1]) %>% ungroup() %>% 
  mutate(index = row_number(),name = paste0(">", name)) %>% 
  filter(len > 0.5*reflen) %>% 
  gather(seq,name,key = "type",value = "fasta") %>% 
  arrange(index, type) %>% select(key,fasta)
tes4 %>% select(fasta) %>% write_tsv("/media/pericles/TEfind/newTE/selectTE.fasta",col_names = FALSE)

tes3 <- bind_rows(TElist %>%  mutate(name2 = paste0(key,"%","embl","%",row_number()),type = "ensemblTE"),
                  pbsv_TE %>%  mutate(name2 = paste0(key,"%","pbsv","%",row_number()),type = "pbsvTE"),
                  rmblast3 %>% mutate(name2 = paste0(key,"%","rmblast","%",row_number()),type = "rmblastTE")) %>% 
  mutate(seq = seq %>% str_to_upper())
tes3 %>% write_tsv("/media/pericles/TEfind/newTE/allTE.tsv")

tes6 <- bind_rows(TElist %>%  mutate(name2 = paste0(key,"%","embl","%",row_number()),type = "ensemblTE"),
                  pbsv_TE %>%  mutate(name2 = paste0(key,"%","pbsv","%",row_number()),type = "pbsvTE"),
                  rmblast3 %>% mutate(name2 = paste0(key,"%","rmblast","%",row_number()),type = "rmblastTE")) %>% 
  mutate(seq = seq %>% str_to_upper(), len = str_length(seq)) %>% 
  group_by(key) %>% mutate(reflen = str_length(seq) %>% .[1]) %>% ungroup() %>% 
  filter(len > 0.5*reflen) 
tes6 %>% write_tsv("/media/pericles/TEfind/newTE/selectTE.tsv")

TElist %>% mutate(seq = seq %>% str_to_upper()) %>% 
  mutate(index = row_number(),name = paste0(">", name)) %>% 
  gather(seq,name,key = "type",value = "fasta") %>% 
  arrange(index, type) %>% select(fasta) %>%
  write_tsv("/media/pericles/TEfind/newTE/refTE.fasta",col_names = FALSE)
TElist %>%  mutate(name2 = paste0(key,"%","embl","%",row_number()),type = "ensemblTE") %>% 
  mutate(seq = seq %>% str_to_upper()) %>% 
  write_tsv("/media/pericles/TEfind/newTE/refTE.tsv")

# pbsv bed generation
pbsv_repmask2 <- read_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_mapped.tsv")
pbsv_bed2 <- pbsv_out %>% 
  left_join(pbsv_ins, by= "name2") %>% 
  mutate(name3 = paste0(chr, "|",start,"|",end,"|",strand,"|pbsv"),len = str_length(seq)) %>%
  group_by(name2,name3) %>% arrange(desc(len)) %>% mutate(fil = row_number()) %>% ungroup() %>% 
  filter(fil ==1) %>% rename(key = TEname) %>% 
  mutate(name3 = paste0(key,"%","pbsv","%",row_number())) %>% 
  select(name2,name3) %>% 
  left_join(pbsv_repmask2,by = "name2")
pbsv_bed2 %>% write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_mapped2.tsv")
pbsv_bed2 %>% mutate(name3 = name3 %>% str_replace("-","\\|") %>% str_remove("-")) %>% 
  select(name3,insstart,insend, TEname, score, strand) %>% 
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_mapped2.bed", col_names = FALSE)

# 
pbsv_anno <- pbsv_out %>% 
  left_join(pbsv_ins, by= "name2") %>% 
  mutate(name3 = paste0(chr, "|",start,"|",end,"|",strand,"|pbsv"),len = str_length(seq)) %>%
  group_by(name3) %>% arrange(desc(len)) %>% mutate(fil = row_number()) %>% ungroup() %>% 
  filter(fil ==1) %>% rename(key = TEname, name = name3) %>% 
  mutate(name3 = paste0(key,"%","pbsv","%",row_number())) %>% 
  select(name,name2,name3) 
pbsv_anno %>% write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_ins_name.tsv")
  
  

