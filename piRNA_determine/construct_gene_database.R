# database make
# /usr/bin/Rscript --silent --slave --vanilla /media/pericles/CLASH/database/construct_gene_database.R;
#######################################################################################################
# gene_TE
#######################################################################################################
# setwd("/media/pericles/CLASH/database/")
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
setwd("/media/pericles/CLASH/database/")
suppressMessages(suppressWarnings(require(tidyverse)))

miRNA <- read_tsv("reffasta/miRNA.tsv", col_names = c("name", "seq")) %>% 
  mutate(type = "miRNA") %>% 
  select(type, name, seq) %>% mutate(name2 = name) %>% 
  separate(name2, sep = ";", into = c("namein", "loc", "ID", "symbol")) %>% 
  separate(namein, sep=" ", into = c("name3", "type")) %>% 
  mutate(type = str_replace(type, "type=", ""),
         type = str_replace(type, ";", "") %>% str_replace("pre_miRNA", "premiRNA"),
         symbol= symbol %>% str_replace(" name=", "") %>% str_replace_all("-", "|"))
tx4 <- miRNA %>% filter(type =="premiRNA") %>% 
  mutate(name3  =  paste0(name3, "_", symbol, "_", type)) %>% 
  select(name3, seq)

embl <- read_tsv("/media/pericles/TEfind/newTE/allTE.tsv") %>% 
  mutate(name = name2 %>% str_replace("-", "|") %>% str_replace("_", "|"),
         name3 = paste0(name, "_", name,  "_",type)) %>% 
  select(name3, seq)

ref <- read_tsv("/media/pericles/TEfind/newTE/refTE.tsv") %>% 
  mutate(name = name2 %>% str_replace("-", "|") %>% str_replace("_", "|"),
         name3 = paste0(name, "_", name,  "_",type)) %>% 
  select(name3, seq)
# edit!
gene1 <- read_tsv("reffasta/gene.tsv", col_names = c("name", "seq")) %>% 
  mutate(type ="gene") %>% 
  select(type, name, seq) %>% mutate(name2 = name) %>% 
  separate(name2, sep = ";", into = c("namein", "loc", "ID", "symbol")) %>% 
  separate(namein, sep=" ", into = c("name3", "type")) %>% 
  mutate(type = str_replace(type, "type=", ""),
         type = str_replace(type, ";", ""),
         symbol= symbol %>% str_replace(" name=", "") %>% str_replace_all("-", "|")) 
flamname <-tribble(
  ~name3, ~symbol, 
  "FBgn9999996", "lncRNA:flamhap1",
  "FBgn9999998", "lncRNA:20Ahap1", 
  "FBgn9999997", "lncRNA:flamhap2", 
  "FBgn9999999", "lncRNA:20Ahap2")

gene1 %>% select(name3, symbol) %>% 
  bind_rows(flamname) %>% 
  write_tsv("anno/dm6_gene_annotation.tsv")
pseudogene1 <- read_tsv("reffasta/pseudogene.tsv", col_names = c("name", "seq")) %>% 
  mutate(type ="pseudogene") %>% 
  select(type, name, seq) %>% mutate(name2 = name) %>% 
  separate(name2, sep = ";", into = c("namein", "loc", "ID", "symbol")) %>% 
  separate(namein, sep=" ", into = c("name3", "type")) %>% 
  mutate(type = str_replace(type, "type=", ""),
         type = str_replace(type, ";", ""),
         symbol = symbol %>% str_replace(" name=", "") %>% str_replace_all("-", "|"))

transcript3 <- bind_rows(gene1 %>% filter(!(str_detect(symbol,"^pre\\|")))) %>%
  mutate(name3  = paste0(name3, "_", symbol, "_", type)) %>% 
  select(name3, seq) 

flam_hap1 <- read_tsv("/media/pericles/TEfind/flam/new_flam/flam_hap1_region.tsv", col_names = c("name", "seq")) %>% 
  mutate(name3 = "FBgn9999996_lncRNA:flamhap1_gene") %>% select(name3,seq)
l20A_hap1 <- read_tsv("/media/pericles/TEfind/flam/new_flam/20A_hap1_region.tsv", col_names = c("name", "seq")) %>% 
  mutate(name3 = "FBgn9999998_lncRNA:20Ahap1_gene") %>% select(name3,seq)
flam_hap2 <- read_tsv("/media/pericles/TEfind/flam/new_flam/flam_hap2_region.tsv", col_names = c("name", "seq")) %>% 
  mutate(name3 = "FBgn9999997_lncRNA:flamhap2_gene") %>% select(name3,seq)
l20A_hap2 <- read_tsv("/media/pericles/TEfind/flam/new_flam/20A_hap2_region.tsv", col_names = c("name", "seq")) %>% 
  mutate(name3 = "FBgn9999999_lncRNA:20Ahap2_gene") %>% select(name3,seq)

hyb_data <- bind_rows(embl, transcript3 %>% filter(!(str_detect(name3, "FBgn0267704|FBgn0287603"))),
                      flam_hap1,l20A_hap1, flam_hap2,l20A_hap2 ) %>% 
  distinct() %>% 
  mutate(name3 = paste0(">", name3),
         index = row_number()) %>% 
  gather(seq, name3, key = "type", value = "fasta") %>% 
  arrange(index, type) %>% select(fasta)
hyb_data %>% write_tsv("dm6_gene.fasta", col_names = FALSE)  
bind_rows(gene1 %>% filter(!(str_detect(symbol,"^pre\\|")))) %>%
  mutate(name3  = paste0(name3, "_", symbol, "_", type)) %>% 
  select(-seq) %>%  
  distinct() %>% write_tsv("anno/dm6_gene.tsv")

hyb_data2 <- bind_rows(ref, transcript3 ) %>% 
  distinct() %>% 
  mutate(name3 = paste0(">", name3),
         index = row_number()) %>% 
  gather(seq, name3, key = "type", value = "fasta") %>% 
  arrange(index, type) %>% select(fasta)
hyb_data2 %>% write_tsv("dm6_gene2.fasta", col_names = FALSE)  

#######################################################################################################
# gene and repeat site
#######################################################################################################
gene2 <- gene1 %>% select(name3, loc) %>% 
  separate(loc, sep = ":", into = c("titania", "miranda")) %>% 
  separate(miranda, sep = "\\.\\.", into = c("start", "end")) %>% 
  mutate(chr = paste0("chr", titania %>% str_remove(" loc=")),
         strand = if_else(str_detect(start, "complement"), "-", "+"),
         start = start %>% str_remove("complement") %>% str_remove("\\(") %>% as.integer(),
         end = end %>% str_remove("\\)"), 
         score = 100) %>% 
  mutate(start = start-1) %>% 
  select(chr, start, end, name3, score, strand) %>% 
  mutate(chr = if_else(chr == "chrmitochondrion_genome", "chrM", chr))
gene2 %>% write_tsv("anno/gene.bed", 
                    col_names = FALSE)



#######################################################################################################
# rRNA & snoRNA ...
#######################################################################################################
rmstRNA <- gene1 %>% 
  filter(str_detect(symbol, "rRNA|tRNA|snRNA|snoRNA|mt:"))

rmstRNA_data <- rmstRNA %>% distinct() %>% 
  mutate(name3 = paste0(">", name3),
         index = row_number()) %>% 
  gather(seq, name3, key = "type", value = "fasta") %>% 
  arrange(index, type) %>% select(fasta)
rmstRNA_data %>% write_tsv("dm6_rmstRNA.fasta", col_names = FALSE)  



########################################################################################################
# annotation name
########################################################################################################
miRNA_name <- read_tsv("reffasta/miRNA.tsv", col_names = c("name", "seq")) %>% 
  mutate(type = "miRNA") %>% 
  select(type, name) %>% mutate(name2 = name) %>% 
  separate(name2, sep = ";", into = c("namein", "loc", "ID", "symbol")) %>% 
  separate(namein, sep=" ", into = c("name3", "type")) %>% 
  mutate(type = str_replace(type, "type=", ""),
         type = str_replace(type, ";", "") %>% str_replace("pre_miRNA", "premiRNA"),
         symbol= symbol %>% str_replace(" name=", "")) %>% filter(type =="premiRNA") %>% 
  mutate(symbol2= symbol %>% str_replace_all("-", "|"),
         name4  =  paste0(name3, "_", symbol2, "_", type)) %>% select(name4,name3,symbol,type)
embl_name <- read_tsv("/media/pericles/TEfind/newTE/allTE.tsv") %>% 
  mutate(name3 = name2 %>% str_replace("-", "|") %>% str_replace("_", "|"),
         name4 = paste0(name3, "_", name3,  "_",type)) %>% separate(name2, sep = "%", into = c("symbol","dataname","index")) %>% 
  select(name4,name3,symbol,type)
gene1_name <- read_tsv("reffasta/gene.tsv", col_names = c("name", "seq")) %>% 
  mutate(type ="gene") %>% 
  select(type, name) %>% mutate(name2 = name) %>% 
  separate(name2, sep = ";", into = c("namein", "loc", "ID", "symbol")) %>% 
  separate(namein, sep=" ", into = c("name3", "type")) %>% 
  mutate(type = str_replace(type, "type=", ""),
         type = str_replace(type, ";", ""),
         symbol= symbol %>% str_replace(" name=", "")) %>%
  mutate(symbol2= symbol %>% str_replace_all("-", "|"),
         name4  = paste0(name3, "_", symbol2, "_", type)) %>% 
  select(name4,name3,symbol,type)
pseudogene1_name <- read_tsv("reffasta/pseudogene.tsv", col_names = c("name", "seq")) %>% 
  mutate(type ="pseudogene") %>% 
  select(type, name) %>% mutate(name2 = name) %>% 
  separate(name2, sep = ";", into = c("namein", "loc", "ID", "symbol")) %>% 
  separate(namein, sep=" ", into = c("name3", "type")) %>% 
  mutate(type = str_replace(type, "type=", ""),
         type = str_replace(type, ";", ""),
         symbol= symbol %>% str_replace(" name=", "")) %>%
  mutate(symbol2= symbol %>% str_replace_all("-", "|"),
         name4  = paste0(name3, "_", symbol2, "_", type)) %>% 
  select(name4,name3,symbol,type)

databasename <- bind_rows(miRNA_name,embl_name,gene1_name,pseudogene1_name) %>% 
  write_tsv("dm6_gene_to_symbol.tsv") 

