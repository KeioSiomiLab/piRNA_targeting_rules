
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(gt)))
ppi <- 600
ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}

outdir <- "/home/qqprosperodd/Desktop/Piwi_article_analysis/fig5/"

#######################################################################################################
# MNase-ChIP analysis
#######################################################################################################

sample_order1 <- c("ni","H3K9me3")
sample_order2 <- c("sense","antisense")

total_read <- read_csv2("/media/hermione/MNaseChIP/read_info.csv") %>% 
  mutate(sample = file %>% str_remove("_1.fastq.gz") %>% str_remove("map\\/"),
         num_seqs = num_seqs %>% str_remove_all(",") %>% as.double()) %>% 
  select(-file)

readdir <- "/media/hermione/MNaseChIP/bigwig/"
file_name <- dir(readdir, pattern = "_luc_1.bedgraph$", full.names = TRUE) %>% 
  str_subset("input|ni|H3K9me3") %>% str_subset("sense|anti_")
file_name2 <- file_name  %>%
  str_replace_all("_luc_1.bedgraph$", "") %>%
  str_replace(paste0(readdir, "/"), "")
df <- vector("list", length(file_name))
for (i in seq_along(file_name)) {
  df[[i]] <- read_tsv(file_name[[i]],
                      col_names = c("chr","start","end","cov"))
  df[[i]][["sample"]] <- file_name2[[i]]
}

luc_bigwig <- bind_rows(df) %>%   left_join(total_read,by = "sample") %>% 
  mutate(cov = 1000000*cov/num_seqs) %>% 
  rowwise() %>%
  mutate(site = list(seq(from = start, to = end-1))) %>% unnest_legacy() %>% 
  select(-start,-end) %>% 
  mutate(site = site-442)

luc_bigwig3 <- luc_bigwig %>% separate(sample, sep = "_", into = c("group","target")) %>% 
  mutate(group = group %>% str_replace("anti","antisense")) %>% 
  mutate(group = group %>% factor(levels = sample_order2)) %>% 
  select(-num_seqs) %>% 
  pivot_wider(names_from = "target", values_from="cov") %>% 
  filter(site <= 925 ) %>% 
  mutate(H3K9me3 = H3K9me3/input) %>% 
  select(chr, group,site, H3K9me3) %>% 
  pivot_longer(c(H3K9me3),  names_to = "target", values_to = "cov")

label_name_order <- c("pMT","FLuc","NLuc")
box_anno <-tribble(
  ~chr, ~start, ~end, ~labelname,
  "FLuc",-442,0,"pMT",
  "NLuc",-442,0,"pMT",
  "FLuc",0,925,"FLuc",
  "NLuc",0,925,"NLuc") %>% 
  mutate(labelname = labelname %>% factor(levels = label_name_order))

luc_bigwig3 %>% mutate(le = paste0(group,"_",target)) %>% 
  ggplot(aes(x = site, y = cov, color = group)) +
  geom_rect(data = box_anno, aes(xmin = start, xmax = end,ymin =-0.05, ymax = -0.002, fill = labelname), inherit.aes = FALSE)+
  geom_text(data = box_anno, aes(x = (start+end)/2, y = -0.025, label = labelname), color = "white", inherit.aes = FALSE)+
  geom_step()+ theme_minimal()+
  facet_wrap(~chr, nrow = 1) +scale_color_manual(values = c("dodgerblue2", "orangered"))+
  scale_x_continuous(breaks = c(-400,-200,0,200,400,600,800)) +
  labs(x = "site", y="ChIP(CPM)/input(CPM)",color = "")+ guides(fill = "none")+
  theme(legend.position = c(0.85,0.85), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c('royalblue3', 'orange2','maroon3'))+
  ggsave(paste0(outdir,"sup_fig5_luc_per_input.png"), width =6, height =3, dpi = 300)+
  ggsave(paste0(outdir,"sup_fig5_luc_per_input.pdf"), width =6, height =3)


# H3K4me3
sample_order1 <- c("ni","H3K4me3")
sample_order2 <- c("sense","antisense")

file_name3 <- dir(readdir, pattern = "_luc_1.bedgraph$", full.names = TRUE) %>% 
  str_subset("input|ni|H3K4me3") %>% str_subset("sense|anti_")
file_name4 <- file_name3  %>%
  str_replace_all("_luc_1.bedgraph$", "") %>%
  str_replace(paste0(readdir, "/"), "")
df2 <- vector("list", length(file_name3))
for (i in seq_along(file_name3)) {
  df2[[i]] <- read_tsv(file_name3[[i]],
                      col_names = c("chr","start","end","cov"))
  df2[[i]][["sample"]] <- file_name4[[i]]
}

luc_bigwig4 <- bind_rows(df2) %>%   left_join(total_read,by = "sample") %>% 
  mutate(cov = 1000000*cov/num_seqs) %>% 
  rowwise() %>%
  mutate(site = list(seq(from = start, to = end-1))) %>% unnest_legacy() %>% 
  select(-start,-end) %>% 
  mutate(site = site-442)

luc_bigwig6 <- luc_bigwig4 %>% separate(sample, sep = "_", into = c("group","target")) %>% 
  mutate(group = group %>% str_replace("anti","antisense")) %>% 
  mutate(group = group %>% factor(levels = sample_order2)) %>% 
  select(-num_seqs) %>% 
  pivot_wider(names_from = "target", values_from="cov") %>% 
  filter(site <= 925 ) %>% 
  mutate(H3K4me3 = H3K4me3/input) %>% 
  select(chr, group,site, H3K4me3) %>% 
  pivot_longer(c(H3K4me3),  names_to = "target", values_to = "cov")

label_name_order <- c("pMT","FLuc","NLuc")
box_anno <-tribble(
  ~chr, ~start, ~end, ~labelname,
  "FLuc",-442,0,"pMT",
  "NLuc",-442,0,"pMT",
  "FLuc",0,925,"FLuc",
  "NLuc",0,925,"NLuc") %>% 
  mutate(labelname = labelname %>% factor(levels = label_name_order))

luc_bigwig6 %>% mutate(le = paste0(group,"_",target)) %>% 
  ggplot(aes(x = site, y = cov, color = group)) +
  geom_rect(data = box_anno, aes(xmin = start, xmax = end,ymin =-0.05, ymax = -0.002, fill = labelname), inherit.aes = FALSE)+
  geom_text(data = box_anno, aes(x = (start+end)/2, y = -0.025, label = labelname), color = "white", inherit.aes = FALSE)+
  geom_step()+ theme_minimal()+
  facet_wrap(~chr, nrow = 1) +scale_color_manual(values = c("dodgerblue2", "orangered"))+
  scale_x_continuous(breaks = c(-400,-200,0,200,400,600,800)) +
  labs(x = "site", y="ChIP(CPM)/input(CPM)",color = "")+ guides(fill = "none")+
  theme(legend.position = c(0.85,0.85), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c('royalblue3', 'orange2','maroon3'))+
  ggsave(paste0(outdir,"sup_fig5_luc_per_input_H3K4me3.png"), width =6, height =3, dpi = 300)+
  ggsave(paste0(outdir,"sup_fig5_luc_per_input_H3K4me3.pdf"), width =6, height =3)


# piRNA expression in vector

piRNA_vector <- read_tsv(paste0("/media/hermione/piRNA_luc/ctrl/vienna/","All_model_vienna.tsv"))



type_order <- c("perfect match","near-perfect match","imperfect match")

piRNA_vector %>% filter(targetname == "NLuc_anti") %>% arrange(sstart2,send2) %>% 
  mutate(cumcount_en = cumsum(CPM), cumcount_st = cumcount_en- CPM) %>% 
  mutate(type = case_when(type == "functional" ~ "near-perfect match",
                          type == "not functional" ~ "imperfect match",
                          TRUE ~ type) %>% factor(levels = type_order)) %>% 
  ggplot(aes(xmin = sstart2-0.5,xmax = send2+0.5,ymin = cumcount_st,ymax = cumcount_en, fill = type))+
  geom_rect() + theme_minimal()+ labs(y = "CPM")+
  annotate("rect", xmin = 0,xmax = 442,ymin =-600, ymax = -0.002, fill = "dodgerblue2")+
  annotate("text", x = 221, y = -300, label = "pMT", color = "white")+
  annotate("rect", xmin = 470,xmax = 1367,ymin =-600, ymax = -0.002, fill = "orangered")+
  annotate("text", x = 919, y = -300, label = "NLuc", color = "white")+
  scale_fill_manual(values = c("maroon1","royalblue1","gray65"))+ guides(fill= "none")+
  theme(axis.title.x = element_blank(),legend.position = "bottom", panel.grid.minor = element_blank())+
  ggsave(paste0(outdir,"piRNA_plasmid_backbone","_cumsum.png"), width =3.5, height =2.3, dpi = 300 )+
  ggsave(paste0(outdir,"piRNA_plasmid_backbone","_cumsum.pdf"), width =3.5, height =2.3)


