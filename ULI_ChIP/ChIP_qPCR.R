#######################################################################################################
# ChIP-qPCR analysis
#######################################################################################################

library(tidyverse)
library(readxl)

mdg1_ChIP <- read_excel("~/ChIP_qPCR/mdg1_repeat_qPCR.xlsx")
ChIP_input <- mdg1_ChIP %>% filter(antibody == "input") %>% 
  rename(Ct2 = Ct) %>% select(-antibody)
mdg1_ChIP2 <- mdg1_ChIP %>% 
  left_join(ChIP_input, by = c("sample","rep","target")) %>% 
  filter(antibody != "input") %>% 
  mutate(per_inp = 5 * (2 ^ (Ct2-Ct)))


sampleorder <- c("sense","anti")
qPCR_make <- function(dataname, outname, antibodyname, fig23_width, height) {
  datafile1 <- dataname
  datafile2 <- datafile1 %>% 
    filter(antibody == antibodyname) %>% 
    mutate(sample = sample %>% factor(levels = sampleorder))
  
  datafile3 <- datafile2 %>% 
    left_join(data_mean %>% filter(sample == "sense") %>% select(-sample), by = c("target")) %>% 
    mutate(scaled_Ct = per_inp/mean_value) 
  
  data_mean2 <- datafile3 %>% group_by(sample,target) %>% 
    mutate(mean_value2 = mean(scaled_Ct),sd=sd(scaled_Ct)) %>% ungroup() %>% 
    select(sample,target,mean_value2,sd) %>% distinct()
  
  g4 <- datafile3 %>% 
    ggplot(mapping = aes(x =sample)) +
    geom_col(data_mean2,mapping = aes(y=mean_value2,fill = sample), color = "black") +
    geom_errorbar(data_mean2,mapping = aes(ymin=mean_value2-sd, ymax=mean_value2+sd), width=.5,
                  position=position_dodge(.9), color = "gray20") + 
    geom_point(mapping = aes(y=scaled_Ct,color = "black"), alpha = 0.9,
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, jitter.height = 0),size = 1.5) +
    theme_minimal()+ guides(color = "none", fill = "none")+ 
    facet_wrap(~ target, scales = "free_y", nrow = 1)+
    scale_fill_manual(values = c("dodgerblue2","orangered"))+
    scale_color_manual(values = c("black"))+
    labs(title = paste0("adjvalue(",outname,")"),
         y = "normalised ChIP score (AU)")+
    #scale_y_continuous(breaks = c(seq(from = 0, to = 1.5,by = 0.1))) +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(size=6))
  g4
  ggsave(paste0("~/ChIP_qPCR/",outname,"_rel2.png"), width =fig23_width, height =height, dpi = 300)
  ggsave(paste0("~/ChIP_qPCR/",outname,"_rel2.pdf"), width =fig23_width, height =height)
}

qPCR_make(dataname = mdg1_ChIP2 %>% filter(target!="Mi-2"), outname = "mdg1_H3K4me3",
          antibodyname = "H3K4me3",fig23_width=4,height=3)
qPCR_make(dataname = mdg1_ChIP2 %>% filter(target!="Mi-2"), outname = "mdg1_H3K9me3",
          antibodyname = "H3K9me3",fig23_width=4,height=3)

