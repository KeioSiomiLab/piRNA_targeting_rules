# database for TEfind
# embl file process
# /usr/bin/Rscript --silent --slave --vanilla TEembl_construct.R database/transposon_sequence_set.embl.txt;
# embldir <- "/media/pericles/TEfind/database/transposon_sequence_set.embl.txt"
#####################################
# embl
####################################
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)
#first TE annotation
embldir <- commandArgs(trailingOnly = TRUE)[1]
embldir <- "/media/pericles/TEfind/database/transposon_sequence_set.embl.txt"
TEannotater <- read_tsv(embldir, col_names = c("V1"))
TEannotater2 <- TEannotater %>% 
  filter((str_detect(V1, "SO_feature") & !str_detect(V1, "RR") & 
            !str_detect(V1, "non_LTR")) | str_detect(V1, "FBte")) %>% 
  mutate(index = if_else(str_detect(V1, "FBte"), 1,0),
         key = cumsum(index)) %>% group_by(key) %>% nest_legacy() %>% 
  mutate(TE = map(data, function(x) {rep(x$V1[[1]], length(x$V1))})) %>% unnest_legacy() %>%
  filter(V1 != TE) %>% select(TE,V1)

TEannotater3 <- TEannotater2 %>% 
  mutate(TE2 = TE %>% str_split("[:blank:]"),
         TE3 = map(TE2, function(x) {x %>% str_subset("\\\\")})) %>% 
  select(-TE2) %>% unnest_legacy() %>% mutate(TE3 = TE3 %>% stringi::stri_replace_last_regex("\\.", ""))
CDSannotater <- TEannotater %>% filter(str_detect(V1, "CDS") | str_detect(V1, "name")) %>% 
  mutate(index = if_else(str_detect(V1, "CDS"), 1,0),
         key = cumsum(index)) %>% group_by(key) %>% nest_legacy() %>% 
  mutate(TE = map(data, function(x) {rep(x$V1[[1]], length(x$V1))})) %>% unnest_legacy() %>%
  filter(V1 != TE) %>% select(TE,V1) %>% 
  mutate(V2 = V1 %>% str_split("\\\\"),
         V3 = map(V2, function(x) {x %>% head(-1) %>% str_c(collapse = "\\")}),
         transcript = map(V2, function(x) {x %>% tail(n=1)})) %>% 
  select(-V2) %>% unnest_legacy() %>% 
  mutate(V3 = V3 %>% str_replace("FT                   /name=", ""),
         V3 = V3 %>% str_replace_all('"', ''),
         V3 = if_else(!str_detect(V3, "\\\\"), paste0("Dmel\\", V3), V3),
         transcript = transcript %>% str_replace('"', '')) %>% 
  select(-V1) %>% rename(V1= TE, TE3 = V3)
#I used stri_replace_last_regex, because "17.6." should be "17.6"
TEannotater4 <- TEannotater3 %>% left_join(CDSannotater, by = c("V1", "TE3")) %>% 
  mutate(V2 = V1 %>% str_split("[:blank:]"),
         V3 = map(V2, function(x) {x %>% str_subset("SO:")}),
         feature = map(V2, function(x) {x %>% .[10]})) %>% 
  select(-TE, -V1,-V2) %>% unnest_legacy()
TEannotater5 <- TEannotater4 %>% separate_rows(V3, sep = ":") %>% 
  filter(str_detect(V3, "\\.")) %>% 
  separate_rows(V3, sep=",") %>% 
  separate(V3, sep = "\\.+", into = c("start", "end")) %>% 
  mutate(feature = if_else(is.na(transcript), feature, paste(feature, transcript, sep = "_")),
         start = start %>% str_replace_all("[:punct:]", "") %>% str_replace("join|complement", "") %>% str_replace("<", ""),
         end = end %>% str_replace_all("[:punct:]", "") %>% str_replace(">", ""),
         TE3 = TE3 %>% str_replace("\\\\","-")) %>% select(-transcript) %>% 
  filter(start !="" & end !="") %>% 
  filter(!(TE3 =="Dmel-Fw2" & feature == "CDS_ORF2"))
TEannotater5 %>% mutate(start = as.double(start) -1, score = ".", strand = "+") %>% 
  write_tsv("/media/pericles/TEfind/database/annotateTE/TE_bed6_full.bed", col_names = FALSE)
TEannotater5 %>% mutate(start = as.double(start) -1, score = ".", strand = "+") %>% 
  filter(!str_detect(feature, "codon|polyA|TATA|site")) %>% 
  write_tsv("/media/pericles/TEfind/database/annotateTE/TE_bed6.bed", col_names = FALSE)
#synonym detect
synonym <- TEannotater %>% filter(str_detect(V1, "synonym") | str_detect(V1,"FBte")) %>% 
  mutate(index = if_else(str_detect(V1, "FBte"), 1,0),
         key = cumsum(index)) %>% group_by(key) %>% nest_legacy() %>% 
  mutate(TE = map(data, function(x) {rep(x$V1[[1]], length(x$V1))})) %>% unnest_legacy() %>%
  filter(V1 != TE) %>% select(TE,V1) %>% 
  mutate(TE2 = TE %>% str_split("[:blank:]"),
         TE3 = map(TE2, function(x) {x %>% str_subset("\\\\")})) %>% 
  select(-TE2) %>% unnest_legacy() %>% mutate(TE3 = TE3 %>% stringi::stri_replace_last_regex("\\.", "")) %>% 
  mutate(V1 =V1 %>% str_replace("SY   synonym: ", "")) %>% 
  select(TE3, V1)
synonym %>% write_tsv("database/annotateTE/TE_synonym.tsv")

#fasta

fasta <- TEannotater %>% filter(!str_detect(V1, "FT|XX|SQ|ID|SY|^AC  |//|^CC  ") | str_detect(V1,"FBte")) %>% 
  mutate(index = if_else(str_detect(V1, "FBte"), 1,0),
         key = cumsum(index)) %>% group_by(key) %>% nest_legacy() %>% 
  mutate(TE = map(data, function(x) {rep(x$V1[[1]], length(x$V1))})) %>% unnest_legacy() %>%
  filter(V1 != TE) %>% select(TE,V1) %>% 
  group_by(TE) %>% nest_legacy() %>% 
  mutate(fasta = map(data, function(x) {x$V1 %>% str_c(collapse = "")})) %>% 
  select(-data) %>% unnest_legacy() %>% 
  mutate(fasta = fasta %>% str_replace_all("[:blank:]", "") %>% str_replace_all("[:digit:]", "") %>% str_to_upper(),
         TE2 = TE %>% str_split("[:blank:]"),
         TE3 = map(TE2, function(x) {x %>% str_subset("\\\\")})) %>% 
  select(-TE2) %>% unnest_legacy() %>% mutate(TE3 = TE3 %>% stringi::stri_replace_last_regex("\\.", "") %>% str_replace_all("\\\\", "-"))
fasta %>% select(TE3, fasta) %>% write_tsv("database/TE_dm.tsv",col_names = FALSE) %>% 
  mutate(index = row_number()) %>% 
  mutate(annotation2 = paste0(">", TE3)) %>% 
  gather(fasta, annotation2, key = "type", value = "fasta") %>% 
  arrange(index, type) %>% select(fasta) %>% 
  write_tsv("database/TE_dm.fasta", col_names = FALSE) 

## generate each TE file:
tes <- fasta %>% select(TE3, fasta) %>% 
  mutate(index = row_number(), TE4 = TE3) %>% 
  mutate(annotation2 = paste0(">", TE3)) %>% 
  gather(fasta, annotation2, key = "type", value = "fasta") %>% 
  arrange(index, type) %>% select(TE4,fasta)

tx <- tes %>% group_by(TE4) %>% nest_legacy()
tx2 <- map(1:length(tx$TE4), function(x) data.frame(tx$data[[x]]))
names(tx2) <- tx$TE4
data <- map(1:length(tx2), function(x) write_tsv(tx2[[x]], 
                                                 paste0("database/eachTE/", names(tx2[x]), ".fasta"),
                                                 col_names = FALSE))

###################################################
# dfam.embl
###################################################
dfamfile <- "/media/pericles/TEfind/database/Dfam.embl"
dfam <- read_tsv(dfamfile, col_names = c("V1"))

fasta_dfam <- dfam %>% filter(!str_detect(V1, "FT|XX|SQ|ID|SY|^AC  |//|^CC  |^CC$|RN|RA|RT|RL|DR|OC|^FH") | str_detect(V1,"^NM  ") | str_detect(V1, "Type:|SubType:|Species:"))
fasta_dfam2 <- fasta_dfam %>% 
  mutate(index = if_else(str_detect(V1, "^NM  "), 1,0),
         key = cumsum(index)) %>% group_by(key) %>% nest_legacy() %>% 
  mutate(TE = map(data, function(x) {rep(x$V1[[1]], length(x$V1))})) %>% unnest_legacy() %>% 
  filter(V1 != TE) %>% select(TE,V1) 
fasta_anno <- fasta_dfam2 %>% filter(str_detect(V1, "Type:|SubType:|Species:")) %>% 
  mutate(type = if_else(str_detect(V1, " Type:"), "family", 
                        if_else(str_detect(V1, " SubType:"), "subfamily",
                                if_else(str_detect(V1, " Species:"), "species",NULL)))) %>% 
  pivot_wider(names_from = "type", values_from = "V1") 
  
anno2 <- fasta_anno %>% 
  mutate(family = family %>% str_remove("CC        Type: "),
         subfamily = subfamily %>% str_remove("CC        SubType:") %>% str_remove("\\s"),
         species = species %>% str_remove("CC        Species: "))
anno2 %>% mutate(TE = TE %>% str_remove("NM   ")) %>% 
  write_tsv("/media/pericles/TEfind/database/class_process_all_species.tsv")

anno3 <- anno2 %>% filter(str_detect(species, "Drosophila")) %>% 
  mutate(TE = TE %>% str_remove("NM   "),
         namekey = TE %>% str_to_lower() %>% str_remove_all("_i|_ltr|_dm|^dm|_rt|-i|-ltr|_fb|-a") %>%
           str_replace("176", "17.6")) %>% 
  arrange(namekey) %>% 
  select(namekey, family, subfamily) %>% distinct()

anno_sub <-tribble(
  ~namekey, ~family,~subfamily,
  "3s18", "LTR",  "Pao",
  "aurora", "LTR",  "Pao",
  "f", "LINE", "I-Jockey",
  "hb","DNA","Unknown",
  "het","LINE", "I-Jockey",
  "jockey","LINE", "I-Jockey",
  "hopper", "DNA","CMC-Transib",
  "p","DNA","P",
  "springer","LTR","Gypsy",
  "flea","LTR","Gypsy",
  "opus","LTR","Gypsy",
  "gate","LTR","Gypsy",
  "x","LINE", "I-Jockey",
  "stalker","LTR","Gypsy",
  "1360","DNA","P",
  "baggins","LINE","R1-LOA",
  "dm88","LTR","Copia",
  "juan","LINE", "I-Jockey",
  "mcclintock","LTR","Gypsy",
  "bari2","DNA","TcMar-Tc1",
  "porto1","LINE", "I-Jockey",
  "tc12","DNA","TcMar-Tc1",
  "dv","Unknown","Unknown",
  "osvaldo","LTR","Gypsy",
  "gandalf","DNA","Unknown",
  "mariner","DNA","TcMar-Tc1",
  "isfun1","DNA","Unknown",
  "bilbo","LINE","R1-LOA",
  "loa","LINE","R1-LOA",
  "uhu","Unknown","Unknown",
  "penelope","Unknown","Unknown",
  "ulysses","LTR","Unknown",
  "tv1","LTR","Unknown",
  "tel","LTR","Unknown",
  "tram","LTR","Unknown",
  "trim","LINE","R1-LOA",
  "paris","Unknown","Unknown",
  "vege","DNA","MITE",
  "mar","Unknown","Unknown",
  "sgm","DNA","MITE",
  "tc3","DNA","TcMar-Tc1",
  "ine1","RC","Helitron",
  "gem","Unknown","Unknown",
  "u","Unknown","Unknown",
  "xanthias","LTR","Copia",
  "fb","DNA","MULE-NOF",
  "hopper2","DNA","CMC-Transib",
  "helitron","RC","Helitron",
  "bungy","SINE","SINE",
  "spock","Unknown","Unknown",
  "worf","Unknown","Unknown",
  "hmsbeagle2","LTR","Gypsy",
  "q","DNA","Q",
  "but1","DNA","Unknown",
  "but2","DNA","Unknown",
  "but3","DNA","Unknown",
  "but4","DNA","Unknown",
  "but5","DNA","Unknown",
  "but6","DNA","Unknown",
  "isbu2","Unknown","Unknown",
  "isbu3","Unknown","Unknown",
  "newton","Unknown","Unknown",
  "galileo","DNA","MULE-NOF",
  "kepler","DNA","MULE-NOF")

# TE seq is flybase
TEs <- fasta %>% select(TE3) %>% mutate(res = "type",TE3 = TE3 %>% str_replace("-","|")) %>% 
  separate(col =TE3, into = c("dis", "name"), sep = "\\|") %>% 
  mutate(namekey = name %>% str_to_lower() %>% str_remove_all("-element|_m|_o|-a$|-b$|-c$")) %>% 
  mutate(namekey=if_else(str_detect(namekey,"^r1"),"r1", namekey) %>% str_remove("-|_"))

withanno <- TEs %>% left_join(bind_rows(anno3, anno_sub), by = "namekey")
withanno %>% arrange(family,subfamily,name) %>% 
  write_tsv("/media/pericles/TEfind/database/class_process.tsv")
# write fasta file for tldr, tebreak
tldr_fasta <- fasta %>% left_join(withanno %>% select(dis,name,family) %>% unite(dis, name,col = "TE3",sep = "-"), by = "TE3")
tldr_fasta %>% mutate(index = row_number()) %>% 
  mutate(annotation2 = paste0(">", family,":",TE3)) %>% 
  gather(fasta, annotation2, key = "type", value = "fasta") %>% 
  arrange(index, type) %>% select(fasta) %>% 
  write_tsv("/media/pericles/TEfind/database/TE_dm_tldr.fasta", col_names = FALSE) 
# extract seq from dfam
fasta_dfam3 <- fasta_dfam2 %>% 
  filter(!str_detect(V1, "Type:|SubType:|Species:|^DE  |^KW  |^FH|^FT")) %>% 
  mutate(index = if_else(str_detect(V1, "OS"), 1,0),
         key = cumsum(index)) %>% group_by(key) %>% nest_legacy() %>% 
  mutate(OS = map(data, function(x) {rep(x$V1[[1]], length(x$V1))})) %>% unnest_legacy() %>% 
  filter(V1 != OS) %>% select(TE,OS,V1) %>% 
  group_by(TE,OS) %>% nest_legacy() %>% 
  mutate(fasta = map(data, function(x) {x$V1 %>% str_c(collapse = "")})) %>% 
  select(-data) %>% unnest_legacy() %>% 
  mutate(fasta = fasta %>% str_replace_all("[:blank:]", "") %>% str_replace_all("[:digit:]", "") %>% str_to_upper(),
         TE = TE %>% str_remove("NM   "),
         OS = OS %>% str_remove("OS   "))
fasta_dfam3 %>% write_tsv("/media/pericles/TEfind/database/dfam_allTEfasta.tsv")
fasta_dfam3 %>% 
  select(TE,fasta) %>% 
  mutate(index = row_number()) %>% 
  mutate(annotation2 = paste0(">", TE)) %>% 
  gather(fasta, annotation2, key = "type", value = "fasta") %>% 
  arrange(index, type) %>% select(fasta) %>% 
  write_tsv("database/dfam_allTE.fasta", col_names = FALSE) 

LTR <- fasta_dfam3 %>% filter(str_detect(OS,"Drosophila|fly"))

# annotation
CDSdfam <- dfam %>% filter(str_detect(V1, "^FT   CDS|FT                   /product=|^NM   ")) %>% 
  mutate(index = if_else(str_detect(V1, "^NM   "), 1,0),
         key = cumsum(index)) %>% group_by(key) %>% nest_legacy() %>% 
  mutate(TE = map(data, function(x) {rep(x$V1[[1]], length(x$V1))})) %>% unnest_legacy() %>% 
  filter(V1 != TE) %>% select(TE,V1) %>% 
  mutate(index = if_else(str_detect(V1, "^FT   CDS"), "region","cdsname"),
         keyin = if_else(str_detect(V1, "^FT   CDS"), 1,0) %>% cumsum()) %>% 
  pivot_wider(names_from = index,values_from = V1) %>% 
  mutate(region = region %>% str_remove("FT   CDS             "),
         cdsname = cdsname %>% str_remove("FT                   /product=") %>% 
           str_remove_all('"'),
         TE = TE %>% str_remove("NM   ")) %>% 
  separate(region, into = c("start","end"), sep = "\\.\\.")
CDSdfam2 <- CDSdfam %>% 
  mutate(cdsname2 = case_when(str_detect(cdsname, "_gag") ~ "gag",
                              str_detect(cdsname, "_pol") ~ "pol",
                              str_detect(cdsname, "_env") ~ "env",
                              str_detect(cdsname, "_tp") ~ "tp",
                              str_detect(cdsname, "_p1") ~ "p1",
                              str_detect(cdsname, "_p2") ~ "p2",
                              str_detect(cdsname, "_p3") ~ "p3",
                              str_detect(cdsname, "_p4") ~ "p4",
                              str_detect(cdsname, "_1p") ~ "1p",
                              str_detect(cdsname, "_2p") ~ "2p",
                              str_detect(cdsname, "_pro") ~ "pro",
                              str_detect(cdsname, "_yr") ~ "yr",
                              str_detect(cdsname, "_orf3") ~ "orf3",
                              str_detect(cdsname, "_ex") ~ "ex",
                              str_detect(cdsname, "_hel") ~ "hel",
                              str_detect(cdsname, "_pw") ~ "pw",
                              str_detect(cdsname, "_px") ~ "px",
                              str_detect(cdsname, "_py") ~ "py",
                              str_detect(cdsname, "_pz") ~ "pz",
                              TRUE ~ cdsname))
CDSdfam2 %>% select(TE, start,end, cdsname2) %>% 
  write_tsv("/media/pericles/TEfind/database/dfam_cds_region.tsv")

###############################################
# GTF construction
###############################################
fasta <- read_tsv("/media/pericles/TEfind/database/TE_dm.tsv",
                  col_names = c("chr","seq")) 

GTF <- fasta %>% mutate(end = str_length(seq),start = 1,type = "gene",strand = "+",
                        score1 = ".",score2 = ".",source = "embl",
                        name = paste0('gene_id "',chr,'"; gene_symbol "',chr,'"')) %>% 
  select(chr, source,type,start,end,score1,strand,score2,name)
GTF %>% write.table("/media/pericles/CLASH/database/gtf/dm6_TE.gtf",
                    col.names = FALSE,quote=FALSE,sep = "\t", row.names = FALSE)


