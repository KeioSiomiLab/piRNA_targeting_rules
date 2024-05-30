
library(tidyverse)
library(readxl)


# target comparison 

luc_make2 <- function(data, outname,target_factor,control,control2,dataset,
                     fig1_width,fig23_width,height, angle_mode = FALSE, luc_mode = TRUE, sample_mode = TRUE,outtable_mode = FALSE){
  exp_data <- read_excel(paste0("~/meeting/",dataset,"_dataset.xlsx"))
  
  target_factor <- target_factor
  data2 <- data %>% 
    pivot_longer(cols = -c("sample","measure","exp"), names_to = "target",values_to = "value") %>% 
    left_join(exp_data, by = "target") %>% 
    mutate(target = target %>% factor(levels = target_factor))
  
  if(luc_mode){
    data3 <- data2 %>% 
      pivot_wider(names_from = "measure",values_from = "value") %>% 
      mutate(adj_value = NLuc/FLuc)
  } else {
    data3 <- data2 %>% 
      pivot_wider(names_from = "measure",values_from = "value") %>% 
      mutate(adj_value = FLuc/NLuc)
  }
  if (sample_mode) {
    target_mean <- data3 %>% filter(ctrl == control) %>% group_by(group) %>% 
      summarise(meanval = mean(adj_value)) %>% ungroup()
  }else{
    target_mean <- data3 %>% filter(target == control) %>% group_by(group) %>% 
      summarise(meanval = mean(adj_value)) %>% ungroup()
  }
  
  data4 <- data3 %>% left_join(target_mean,by = c("group")) %>% 
    mutate(rel_value = adj_value / meanval) %>% 
    group_by(target, ctrl,group) %>% mutate(mean_value = mean(rel_value),sd = sd(rel_value)) %>% ungroup() 
  
  data4_mean <- data4 %>% select(target,ctrl,group,mean_value,sd) %>% distinct()
  
  g2 <- data4 %>% 
    ggplot(mapping = aes(x =target)) +
    geom_col(data4_mean,mapping = aes(y=mean_value,fill = ctrl), color = "black") +
    geom_errorbar(aes(ymin=mean_value-sd, ymax=mean_value+sd), width=.5,
                  position=position_dodge(.9), color = "gray20") + 
    geom_point(mapping = aes(y=rel_value,color = "black"),
               alpha = 0.9,position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, jitter.height = 0),size = 1.5) +
    theme_minimal()+ guides(color = "none", fill = "none")+
    scale_fill_manual(values = c("dodgerblue2","orangered"))+
    scale_color_manual(values = c("black"))+
    labs(title = paste0("adjvalue(",outname,")"),
         y = "luciferase activity (relative to target)")+
    scale_y_continuous(breaks = c(seq(from = 0, to = 1.3,by = 0.1))) +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(size=6))
  g2
  if(angle_mode){
    g2 <- g2 + theme(axis.text.x = element_text(angle = 30,hjust = 1))
  }
  g2
  ggsave(paste0("~/meeting/",outname,"_rel.png"), width =fig23_width, height =height, dpi = 300)
  ggsave(paste0("~/meeting/",outname,"_rel.pdf"), width =fig23_width, height =height)
  
}

luc_make2(data = data210806_locus, outname = "data210806_locus_dataset",dataset ="210806_locus",
         target_factor = c("5UTR_sense","5UTR_anti","3UTR_sense","3UTR_anti",
                           "Rpl23intron_sense","Rpl23intron_anti","rasintron_sense","rasintron_anti"),
         control = "control",control2 = "target", fig1_width = 5,fig23_width = 3.5, height =3.5, angle_mode = TRUE)

luc_make2(data = data210809_flam, outname = "data210809_flam_dataset",dataset = "210809_flam",
         target_factor = c("flam0.6_F","flam0.6_R","flam0.3_F","flam0.3_R","flam0.2_F","flam0.2_R",
                           "flam0.1_F","flam0.1_R","flam0.08_F","flam0.08_R"),
         control = "control",control2 = "target", fig1_width = 10,fig23_width = 4, height =3.5,
         angle_mode = TRUE,luc_mode = FALSE)

luc_make2(data = data210831_KD, outname = "data210831_KD_dataset",dataset = "210831_KD",
         target_factor = c("siEGFP_sense","siEGFP_anti","siPIWI_sense","siPIWI_anti",
                           "siGTSF1_sense","siGTSF1_anti","siNXF2_sense","siNXF2_anti"),
         control = "control", control2 = "target",
         fig1_width = 5,fig23_width = 3, height =3, angle_mode = TRUE)

luc_make2(data = data210902_other, outname = "data210902_other_dataset",dataset = "210902_other",
         target_factor = c("mdg1_sense","mdg1_anti","piRNA134303_sense",
                           "piRNA134303_anti","piRNA268036_sense","piRNA268036_anti"),
         control = "control", control2 = "target",
         fig1_width = 5,fig23_width = 3, height =3.5,angle_mode = TRUE)

