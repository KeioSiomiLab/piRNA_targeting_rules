# GTF database make
# /usr/bin/Rscript --silent --slave --vanilla /media/pericles/CLASH/GTF_process.R;
#######################################################################################################
# annotation
#######################################################################################################
# setwd("/media/pericles/CLASH/database/anno/")
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

suppressMessages(suppressWarnings(require(tidyverse)))
gtf <- read_tsv("/media/pericles/CLASH/database/gtf/dm6.gtf",
                col_names= c("chr", "sourse", "type", "start", "end", "length", "direction", "dis", "gene_name"),
                col_types = cols(chr = col_character(),
                                 sourse = col_character(),
                                 type = col_character(),
                                 start = col_double(),
                                 end = col_double(),
                                 length = col_character(),
                                 direction = col_character(),
                                 dis = col_character(),
                                 gene_name = col_character()))
gtf2 <- gtf %>% separate(gene_name, into = c("name","genesymbol","transcript","tfsymbol"), sep = "\\;") %>% 
  select(-length, -dis) %>% 
  mutate(name = name %>% str_remove("gene_id") %>% str_remove_all("[:punct:]")%>% str_remove_all("[:space:]"),
         transcript = transcript %>%  str_remove("transcript_id") %>% str_remove_all("[:punct:]")%>% str_remove_all("[:space:]"),
         genesymbol = genesymbol %>%  str_remove("gene_symbol") %>% str_remove_all('"')%>% str_remove_all("[:space:]"),
         tfsymbol = tfsymbol %>%  str_remove("transcript_symbol") %>% str_remove_all('"')%>% str_remove_all("[:space:]"),
         chr = paste0("chr", chr)) 
gtf2 %>% write_tsv("/media/pericles/CLASH/database/gtf/dm6_use.tsv")

mdg4 <- gtf %>% 
  filter(str_detect(gene_name, "mdg4")) %>% 
  mutate(trans = str_extract(gene_name,"FBtr\\d+") %>% str_replace_na(replacement = "gene"),
         gene = str_extract(gene_name,"FBgn\\d+")) %>% 
  group_by(gene,trans) %>% mutate(ariel = paste0(direction %>% unique(),collapse = "")) %>% 
  ungroup() %>% 
  filter(ariel !=".-+" & ariel !="+") %>% 
  select(-trans,-gene,-ariel)

gtf %>%   filter(!(str_detect(gene_name, "mdg4"))) %>% 
  filter(chr %in% c("2L","2R","3L","3R","4","X","Y")) %>% 
  bind_rows(mdg4) %>% 
  mutate(chr = paste0("chr", chr)) %>% 
  write.table("/media/pericles/CLASH/database/gtf/dm6_use.gtf",col.names = FALSE,quote=FALSE,sep = "\t", row.names = FALSE)

tes <- gtf %>% mutate(chr = paste0("chr", chr)) %>% 
  filter(chr %in% paste0("chr", c("2L","2R","3L","3R","4","X","Y"))) %>% 
  mutate(id = gene_name %>% str_extract('gene_symbol "[:graph:]+";') %>% str_replace("gene_symbol","gene_id"),
         tf = gene_name %>% str_extract('transcript_symbol "[:graph:]+";') %>% str_replace("transcript_symbol","transcript_id"),
         name = if_else(type == "gene", id, paste0(id, " ",tf))) %>% select(-id,-tf,-gene_name) %>% 
  write.table("/media/pericles/CLASH/database/gtf/dm6_igv.gtf",col.names = FALSE,quote=FALSE,sep = "\t", row.names = FALSE)
gene <- gtf2 %>% 
  filter(type =="gene") %>% 
  mutate(gene_start = start,gene_end = end) %>%
  select(name, gene_start,gene_end)
exon <- gtf2 %>% 
  filter(type =="exon") %>% filter(name != "FBgn0002781") %>% 
  left_join(gene, by = "name")
exon2 <- exon %>% 
  mutate(new_start = if_else(direction=="+", start - gene_start, gene_end-end),
         new_end = if_else(direction=="+", end - gene_start,gene_end-start)) %>% 
  select(name,transcript, new_start,new_end)
intron <- exon2 %>% 
  group_by(name,transcript) %>% nest_legacy() %>% 
  mutate(start = map(data, function(x) {start = head(x$new_end-1, -1)}),
         end = map(data, function(x) {end = tail(x$new_start-1, -1)})) %>% 
  select(-data) %>% unnest_legacy() 
anno <- read_tsv("/media/pericles/CLASH/database/anno/dm6_gene_annotation.tsv") %>% 
  rename(name = name3)
intron %>% left_join(anno, by = "name") %>% 
  mutate(direction2 = "+", name = paste0(name, "_",symbol, "_gene")) %>% 
  select(name, start, end, direction2) %>% 
  write_tsv("/media/pericles/CLASH/database/anno/dm6_gene_splicesite.txt", col_names = FALSE)

#######################################################################################################
# get gene type
#######################################################################################################
suppressMessages(suppressWarnings(require(tidyverse)))
gtf <- read_tsv("/media/pericles/CLASH/database/gtf/dm6.gtf",
                col_names= c("chr", "sourse", "type", "start", "end", "length", "direction", "dis", "gene_name"),
                col_types = cols(chr = col_character(),
                                 sourse = col_character(),
                                 type = col_character(),
                                 start = col_double(),
                                 end = col_double(),
                                 length = col_character(),
                                 direction = col_character(),
                                 dis = col_character(),
                                 gene_name = col_character()))
gettype <- gtf %>% filter(!(type %in% c("CDS","gene","3UTR","5UTR","start_codon","stop_codon","exon"))) %>% 
  separate(gene_name, into = c("name","genesymbol","transcript","tfsymbol"), sep = "\\; ") %>% 
  select(type, name, genesymbol) %>% distinct() %>% 
  mutate(name = name %>% str_remove('gene_id "') %>% str_remove('"'),
         genesymbol = genesymbol %>% str_remove('gene_symbol "') %>% str_remove('"'))

flamname <-tribble(
  ~type, ~name, ~genesymbol,
  "ncRNA","FBgn9999996", "lncRNA:flamhap1",
  "ncRNA","FBgn9999998", "lncRNA:20Ahap1", 
  "ncRNA","FBgn9999997", "lncRNA:flamhap2", 
  "ncRNA","FBgn9999999", "lncRNA:20Ahap2")
gettype %>% bind_rows(flamname) %>% 
  write_tsv("/media/pericles/CLASH/database/gtf/dm6_RNAtype.tsv")

#######################################################################################################
# TSS site
#######################################################################################################
tss <- gtf2 %>% filter(chr %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY")) %>% 
  filter(!(type %in% c("exon","start_codon","3UTR","5UTR", "CDS","stop_codon","gene"))) %>% 
  select(chr,type, start, end, direction, name,genesymbol,transcript, tfsymbol) %>% 
  mutate(tss = if_else(direction == "+", start-1, end)) %>% 
  write_tsv("/media/pericles/CLASH/database/gtf/dm6_tss.tsv")
tss %>%   mutate(score =0, tss = if_else(direction == "+", start-1, end-1), 
                 tss2 = if_else(direction == "+", start, end)) %>% 
  select(chr, tss,tss2, transcript,score, direction ) %>% 
  arrange(chr, tss,tss2) %>% 
  write_tsv("/media/pericles/CLASH/database/gtf/dm6_tss.bed", col_names = FALSE)

# 5UTR for CAGE
UTR5 <- gtf2 %>% filter(chr %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY")) %>% 
  filter(type == "5UTR") %>% 
  select(chr,type, start, end, direction, name,genesymbol,transcript, tfsymbol) %>% 
  write_tsv("/media/pericles/CLASH/database/gtf/dm6_5UTR.tsv")
UTR5 %>%  mutate(score =0) %>% 
  select(chr, start,end, transcript,score, direction ) %>% 
  arrange(chr, start,end) %>% 
  write_tsv("/media/pericles/CLASH/database/gtf/dm6_5UTR.bed", col_names = FALSE)

