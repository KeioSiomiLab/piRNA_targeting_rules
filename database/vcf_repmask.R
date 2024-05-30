# vcf repmask
# /usr/bin/Rscript --silent --slave --vanilla vcf_repmask.R;
#######################################################################################################
# pbsv
#######################################################################################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gt)))
suppressMessages(suppressWarnings(require(ggpointdensity)))
suppressMessages(suppressWarnings(require(viridis)))

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
TEgroup <- read_tsv("/media/pericles/TEfind/database/class_process.tsv") %>%
  unite(dis, name, sep = "-",col = "name") %>%
  select(name, family, subfamily)

pbsv_bed <- read_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins.bed",
                     col_names = c("chr","start","end","name2","score2","strand2"))

filedir <- "/media/pericles/TEfind/vcfprocess/rmblast/pbsv_vcf_ins_read.fasta.out.gff"
anno_pbsv <-  read_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins_annotation.tsv")

pbsv_repmask <- read_tsv(filedir,  col_names = c("name2","source", "dis1","insstart","insend", "score","strand","dis2","name"),skip = 3)

# filteration of imprecise insertion!

imprecise_locus <- pbsv_repmask %>%
  separate(name,into = c("dis3","TEname","TEstart","TEend"), sep = " ") %>%
  mutate(TEname = TEname %>% str_remove_all('"|Motif:')) %>%
  left_join(pbsv_bed %>% select(chr, start,end,name2), by= "name2") %>% 
  mutate(ID = name2 %>% str_extract("^pbsv.INS.\\d+")) %>% left_join(anno_pbsv, by = "ID") %>% filter(certainty == "imprecise") %>% 
  select(chr, start,end, TEname,score, strand) %>% arrange(chr,start,end) %>% 
  write_tsv("/media/pericles/TEfind/vcfprocess/pbsv_vcf_ins_imprecise_locus.bed", col_names = FALSE)

pbsv_repmask2 <- pbsv_repmask %>%
  separate(name,into = c("dis3","TEname","TEstart","TEend"), sep = " ") %>%
  mutate(TEname = TEname %>% str_remove_all('"|Motif:')) %>%
  left_join(pbsv_bed %>% select(chr, start,end,name2), by= "name2") %>% 
  mutate(ID = name2 %>% str_extract("^pbsv.INS.\\d+")) %>% left_join(anno_pbsv, by = "ID") %>% 
  filter(certainty == "precise")
pbsv_repmask2 %>% write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_mapped_all.tsv")
refTE <- pbsv_repmask2 %>%
  select(name2) %>% distinct() %>%
  mutate(length = name2 %>% str_extract("@\\d+") %>% str_remove("@") %>% as.integer(),
         start = 0)



pbsv_repmask2 %>% select(name2, insstart,insend,score,strand,TEname,TEstart,TEend,chr,start,end) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_mapped.tsv")
pbsv_repmask2 %>% mutate(name3 = paste0(name2,"&", insstart,"&",insend,"%",chr,"%",start,"%",end,"%",strand),
                         TEstart = as.integer(TEstart) -1) %>%
  select(TEname,TEstart,TEend,name3) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_locus.bed",
            col_names = FALSE)


pbsv_repmask3 <- pbsv_repmask2 %>%
  mutate(len = insend-insstart, score2 = score*len) %>%
  group_by(chr, start, end, name2, strand,TEname) %>%
  summarise(score3 = sum(score2), total_len = sum(len)) %>% ungroup() %>%
  mutate(new_score = score3/total_len) %>%
  group_by(chr, start,end) %>% filter(new_score == min(new_score)) %>% filter(total_len == max(total_len)) %>%  ungroup() %>%
  select(chr, start, end, TEname, new_score, strand,total_len) %>% distinct() %>%
  mutate(type = "pbsv")

pbsv_repmask3 %>% select(chr, start, end, TEname, new_score, strand) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins.bed",col_names = FALSE)

pbsv_repmask2 %>%
  mutate(len = insend-insstart, score2 = score*len) %>%
  group_by(chr, start, end, name2, strand,TEname) %>%
  summarise(score3 = sum(score2), total_len = sum(len)) %>% ungroup() %>%
  mutate(new_score = score3/total_len) %>%
  group_by(chr, start,end) %>% filter(new_score == min(new_score)) %>% filter(total_len == max(total_len)) %>%  ungroup() %>%
  mutate(end = start+1) %>%
  select(chr, start, end, name2, new_score, strand) %>% distinct() %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_ref.bed",col_names = FALSE)

pbsv_repmask3 %>% mutate(TEname = paste0(TEname,"@pbsv")) %>%
  select(chr, start, end, TEname, new_score, strand) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_venn.bed",col_names = FALSE)
pbsv_repmask3 %>% select(chr, start, end, TEname, new_score, strand) %>%
  mutate(length = end-start) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins.tsv")

pbsv_out <- pbsv_repmask2 %>%
  mutate(len = insend-insstart, score2 = score*len) %>%
  group_by(chr, start, end, name2, strand,TEname) %>%
  summarise(score3 = sum(score2), total_len = sum(len)) %>% ungroup() %>%
  mutate(new_score = score3/total_len) %>%
  group_by(chr, start,end) %>% filter(new_score == min(new_score)) %>% filter(total_len == max(total_len)) %>%  ungroup() %>%
  select(chr, start, end, TEname, new_score, strand,total_len,name2) %>% distinct() %>%
  mutate(type = "pbsv") %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_ins_anno.tsv")


rmblast_len <- read_tsv("/media/pericles/TEfind/rmblast/dm6_rmblast_mask.tsv") %>%
  rename(TEname = name) %>% mutate(TElen = TEend-TEstart)
rmblast_repmask <- read_tsv("/media/pericles/TEfind/vcfprocess/embl_pbsv_del_repeat.bed",
                            col_names = c("chr","start", "end", "TEname","new_score","strand")) %>%
  mutate(type = "rmblast") %>%
  left_join(rmblast_len %>% select(chr, start,end,strand, TEname,TElen), by = c("chr","start","end","TEname","strand"))

rmblast_mapping <- read_tsv("/media/pericles/TEfind/rmblast/dm6_rmblast_mask.tsv")

rmblast_locus <- rmblast_repmask %>%
  left_join(rmblast_mapping %>% rename(TEname = name, new_score = score),
            by = c("chr","start","end","TEname","new_score","strand"))
rmblast_locus %>%
  mutate(name3= paste0("rmblast","%",chr,"%",start,"%",end,"%",strand),
         TEstart = as.integer(TEstart) -1) %>%
  select(TEname, TEstart,TEend,name3) %>%
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_rmblast_locus.bed",
            col_names = FALSE)

tree_repmask <- bind_rows(pbsv_repmask3 %>% rename(TElen = total_len), rmblast_repmask) %>%
  left_join(TEgroup %>% rename(TEname = name), by = "TEname") %>%
  mutate(TEname = paste0(TEname,"@",type))

gafa2 <- tree_repmask %>% separate(TEname, sep = "@",into =c("TEname","dis")) %>%
  group_by(type, family,TEname) %>% summarise(count = n()) %>% ungroup() %>%
  group_by(type) %>%
  arrange(desc(count)) %>%
  top_n(20) %>%
  gt() %>%
  gtsave("/media/pericles/TEfind/vcfprocess/fig/count.pdf")

TE_len <- read_tsv("/media/pericles/TEfind/database/TE_dm.fasta.fai",
                   col_names = c("TEname","len","dis1","dis2","dis3")) %>%
  select(TEname, len)

tree_repmask2 <-  tree_repmask %>% select(chr, start,end,TEname,new_score,strand,TElen) %>%
  separate(TEname, sep = "@",into =c("TEname","type")) %>%
  left_join(TE_len,by = "TEname") %>%
  mutate(fullness = (TElen)/len,
         anno = case_when(fullness > 1.3 ~ "multi",
                          fullness > 0.9 ~ "full",
                          fullness > 0.6 ~ "almost_full",
                          fullness > 0.4 ~ "short",
                          fullness > 0.2 ~ "little",
                          TRUE ~ "dust"))


tree_repmask2 %>% mutate(name = paste0(TEname,"@",type,"@",anno)) %>%
  select(chr, start,end,name,new_score,strand) %>%
  write_tsv("/media/pericles/TEfind/newTE/OSCrepeat_pbsv.bed",
            col_names = FALSE)


# insertion annotation

homohetero <- pbsv_repmask3 %>% 
  left_join(pbsv_repmask2 %>% select(chr,start, end, TEname,strand,genotype) %>% distinct(), by = c("chr","start", "end", "TEname","strand"))

homohetero %>% 
  filter(genotype =="homo") %>% 
  select(chr, start, end, TEname, new_score, strand,genotype) %>% arrange(chr,start, end) %>% 
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_homo_fil.bed",
            col_names = FALSE)
homohetero %>% 
  filter(genotype =="hetero") %>% 
  select(chr, start, end, TEname, new_score, strand,genotype) %>% arrange(chr,start, end) %>% 
  write_tsv("/media/pericles/TEfind/vcfprocess/rmblast/dm6_pbsv_hetero_fil.bed",
            col_names = FALSE)








