#######################################################################################################
# separate chimera file into piRNA and target sequence
#######################################################################################################
# /usr/bin/Rscript  piclash_make_vienna2.R Process3/{/.}_chimera.tsv vienna/{/.}.rnaplex vienna/{/.}
# hybdir <- "/media/hermione/SMEDWI3/Process3/SRR8842976_chimera.tsv"
# rnaplexname <- "/media/hermione/SMEDWI3/vienna/SRR8842976.rnaplex"
# outname <- "/media/hermione/SMEDWI3/vienna/SRR8842976"


.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")
library(tidyverse)
library(RRNA)
ct2coord_up<-function(input){
  group=0;
  firstrun=1
  ### First NT in the file is at 0,0 ####

  #### Create two fake base-pairs in the begining to make it work ###
  mnp<-min(input$pos)
  map<-max(input$pos)
  a1<-map+1
  a2<-map+2
  b1<-mnp-1
  b2<-mnp-2

  new1<-NULL
  new1$pos[1]=a1
  new1$before[1]=(a1-1)
  new1$after[1]=(a1+1)
  new1$seq[1]="A"
  new1$pos2[1]=a1
  new1$bound[1]=b1

  new2<-NULL
  new2$pos[1]=a2
  new2$before[1]=(a2-1)
  new2$after[1]=(a2+1)
  new2$seq[1]="A"
  new2$pos2[1]=a2
  new2$bound[1]=b2

  new3<-NULL
  new3$pos[1]=b1
  new3$before[1]=(b1-1)
  new3$after[1]=(b1+1)
  new3$seq[1]="A"
  new3$pos2[1]=b1
  new3$bound[1]=a1

  new4<-NULL
  new4$pos[1]=b2
  new4$before[1]=(b2-1)
  new4$after[1]=(b2+1)
  new4$seq[1]="A"
  new4$pos2[1]=b2
  new4$bound[1]=a2

  input<-rbind(new1,new2,input,new3,new4)
  input<-as.data.frame(input)
  mnp<-min(input$pos)

  output<-NULL
  nextNT=input[input$pos==mnp,]
  if(nextNT$bound!=0){
    #### Set Coordinates to (0,0) and (0,1) #####
    output$pos[1]=mnp
    output$x[1]=0
    output$y[1]=0
    output$pos[2]=nextNT$bound
    output$x[2]=0
    output$y[2]=1
  }
  stems<-NULL
  stems[[1]]<-c(output$pos[1],output$pos[2])

  j<-3
  prev<-mnp
  npos<-input$after[input$pos==mnp]
  nextNT=input[input$pos==npos,]
  mp<-max(input$pos)
  while(length(stems)>0){
    newstems<-NULL
    newloops<-NULL
    ns<-1
    nl<-1
    for(i in 1:length(stems)){
      s1<-stems[[i]]
      #### Use stemCoord to generate stem coordinates ###
      p1<-s1[1]
      p2<-s1[2]
      if(firstrun==1){
        x3=-1
        y3=0
        firstrun=0
      }else{
        ### If prev.x < this.x p3=-1 ##
        prev<-input$before[input$pos==p1]
        x3<-output$x[output$pos==prev]
        y3<-output$y[output$pos==prev]
      }
      x1<-output$x[output$pos==p1]
      y1<-output$y[output$pos==p1]
      x2<-output$x[output$pos==p2]
      y2<-output$y[output$pos==p2]
      sout<-stemCords_up(input,p1,p2,x1,y1,x2,y2,x3,y3)
      sdat<-sout[[1]]
      sdat<-as.data.frame(sdat)
      sdat<-sdat[sdat$pos!=p1,]
      sdat<-sdat[sdat$pos!=p2,]
      l<-dim(sdat)[1]
      if(l>0){
        s<-length(output$pos)
        s<-s+1
        for(sd in 1:length(sdat$pos)){
          output$pos[s]<-sdat$pos[sd]
          output$x[s]<-sdat$x[sd]
          output$y[s]<-sdat$y[sd]
          s<-s+1
        }
      }
      newloops[[nl]]=sout[[2]]
      nl<-nl+1
    }
    ##### Loop through all loop starts ###
    newstems<-NULL
    for( i in 1:length(newloops)){
      lps<-newloops[[i]]
      lp1<-loopLength(input,lps[1])
      #### Add stems to newstems ###
      if(length(lp1)>1){
        nstm<-lp1[[2]]
        for(nt in 1:length(nstm)){
          newstems[[ns]]<-nstm[[nt]]
          ns<-ns+1
        }
      }
      #### Get Coordinates for each position and add to output ###
      pp<-input$before[input$pos==lps[1]]
      if(output$x[output$pos==lps[1]] < output$x[output$pos==pp]){
        p3=1
      }else{
        p3=1
      }

      lout<-genCords(lp1,lps[1],lps[2],output,p3)
      pt<-lp1[[1]]
      tout<-lout
      lpl<-length(pt$pos)
      pt$pos=pt$pos[lpl:1]
      lout$pos=pt$pos
      lout<-as.data.frame(lout)
      lout<-lout[lout$pos!=lps[1],]
      lout<-lout[lout$pos!=lps[2],]
      s<-length(output$pos)
      s<-s+1
      for(lo in 1:length(lout$pos)){
        output$pos[s]<-lout$pos[lo]
        output$x[s]<-lout$x[lo]
        output$y[s]<-lout$y[lo]
        output$group[s]=group
        s<-s+1
      }
      group=group+1
    }
    stems<-newstems

  }
  my<-min(output$y)-5
  xy<-max(output$y)+5
  mx<-min(output$x)-5
  xx<-max(output$x)+5
  output<-as.data.frame(output)
  input<-as.data.frame(input)
  all<-merge(output,input,by="pos")
  names(all)[1]="num"
  ### Remove first two and last two ##
  rmax<-dim(all)[1]
  rmax<-rmax-2
  all<-all[3:rmax,]
  return(all)
}
stemCords_up<-function(input,p1,p2,x1,y1,x2,y2,x3,y3){
  ##### INput is a CT file
  ##### P1 is the position of the 5' end of the stem
  ##### P2 is the position of its bound partner
  ##### P3 is the orientation of the stem -1 or 1
  #### Given the first BasePair of a stem it fills in the rest of the stem until
  #### It comes to an unbound member

  output<-NULL
  output$pos[1]=p1
  output$x[1]=x1
  output$y[1]=y1
  output$pos[2]=p2
  output$x[2]=x2
  output$y[2]=y2
  ends<- -1
  #### Determine perpendicular vector ###
  vperp<-NULL
  vperp$x=x3-x1
  vperp$y=y3-y1
  v<-rotateS((x2-x1),(y2-y1),0,0,(pi/2))
  m<-sqrt(v[1]^2+v[2]^2)+sqrt(vperp$x^2+vperp$y^2)
  ang<-acos((v[1]*vperp$x+v[2]*vperp$y)/m)
  if(ang<pi/2){
    v<-rotateS((x2-x1),(y2-y1),0,0,(-1*pi/2))
  }
  vect<-NULL
  vect$x=v[1]
  vect$y=v[2]
  nv<-sqrt((vect$x)^2+(vect$y)^2)
  vect$x<-vect$x/nv
  vect$y<-vect$y/nv
  prev1=p1
  j<-3
  mf<-length(input$pos)+5
  mp<-max(input$pos)
  while(j<mf){
    npos=input$after[input$pos==prev1]
    if(npos<(mp+1)){
      nextNT=input[input$pos==npos,]
      if(dim(nextNT)[1] < 1){
        j=mf+5
      }else{
        if(nextNT$bound==0){
          ##### Position is not bound done ####
          ## Return the last position as Ends #
          bound<-input$bound[input$pos==prev1]
          ends<-c(prev1,bound)
          j=mf+5
        }else{

          pbpos=input$bound[input$pos==prev1]
          pbafter=input$before[input$pos==pbpos]
          bpos=input$bound[input$pos==npos]
          if(pbafter==bpos){
            ### Check that prev bound is on the same chain ##
            output$pos[j]=npos
            output$x[j]=output$x[output$pos==prev1]+vect$x[1]
            output$y[j]=output$y[output$pos==prev1]+vect$y[1]
            j<-j+1
            if(j>3){
              pbpos=input$bound[input$pos==prev1]
              bpos=input$bound[input$pos==npos]
              output$pos[j]=bpos
              output$x[j]=output$x[output$pos==pbpos]+vect$x[1]
              output$y[j]=output$y[output$pos==pbpos]+vect$y[1]
              j<-j+1
            }
          }else{
            ends<-c(prev1,pbpos)
            j=mf+5

          }
          prev1=npos
          #################################################
        }
      }
    }
    else{
      j=mf+5
    }
    ### Second NT is at 0,1 ###
  }
  l<-NULL
  l[[1]]<-output
  l[[2]]<-ends
  return(l)
}

hybdir <- commandArgs(trailingOnly = TRUE)[1]
rnaplexname <- commandArgs(trailingOnly = TRUE)[2]
outname <- commandArgs(trailingOnly = TRUE)[3]

hybtest <- read_tsv(hybdir)
rnaplex <- read_tsv(rnaplexname,
                    col_names = "data") %>%
  mutate(index = case_when(str_detect(data,"piRNA") ~ "piRNA",
                           str_detect(data,">\\d+_\\d+") ~ "readname",
                           TRUE ~ "match"),
         counter = if_else(index =="piRNA", 1,0) %>% cumsum()) %>%
  pivot_wider(names_from = "index",values_from = data) %>%
  drop_na()


rnaplex2 <- rnaplex %>% mutate(match = match %>% str_replace(" \\(", "  \\(") %>% str_replace_all("   ", "  ") %>% str_replace_all("\\) 1", "\\)  1") %>%
                                 str_replace_all("\\. 1", "\\.  1") %>% str_replace_all(" :  ", "  :  ") %>% str_remove(":  ")) %>%
  separate(match, sep = "\\s\\s", into = c("bind", "p3_inno", "p5_inno", "dG_raw")) %>%
  separate(bind, sep = "&", into = c("p3bind", "p5bind")) %>%
  separate(p3_inno, sep = ",", into = c("p3_start","p3_end")) %>%
  separate(p5_inno, sep = ",", into = c("p5_start","p5_end")) %>%
  mutate(p3_start = p3_start %>% str_remove_all(" ") %>% as.integer(),
         p3_end = p3_end %>% str_remove_all(" ") %>% as.integer(),
         p5_start = p5_start %>% str_remove_all(" ") %>% as.integer(),
         p5_end = p5_end %>% str_remove_all(" ") %>% as.integer()) %>%
  mutate(piRNA = piRNA %>% str_remove(">piRNA"), readname = readname %>%  str_remove(">"),
         dG = dG_raw %>% str_extract("\\(.+\\)") %>%  str_remove_all("\\(|\\)") %>% as.double())
rnaplex3 <- hybtest %>%
  left_join(rnaplex2,by = c("readname","piRNA")) %>%
  mutate(piRNAmatch = paste0(strrep(".", p3_start-1), p3bind, strrep(".", pilen - p3_end)) %>% str_replace_all("\\)", "\\("),
         targetmatch = paste0(strrep(".", p5_start-1), p5bind, strrep(".", targetlen - p5_end)) %>% str_replace_all("\\(", "\\)")) %>%
  select(-contains("p3"), -contains("p5"), -contains("dis"), -dG_raw) %>% rename(index = counter) %>%
  mutate(piRNAmatch = if_else(!(str_detect(piRNAmatch, "\\(")), strrep(".", pilen),piRNAmatch),
         targetmatch = if_else(!(str_detect(piRNAmatch, "\\(")), strrep(".", targetlen),targetmatch))
rnaplex3 %>% write_tsv(paste0(outname,"_vienna.tsv"))
#vienna
vienna_struct <- hybtest %>%
  left_join(rnaplex2,by = c("readname","piRNA")) %>%
  mutate(piRNAmatch = paste0(strrep(".", p3_start-1), p3bind, strrep(".", pilen - p3_end)) %>% str_replace_all("\\)", "\\("),
           targetmatch = paste0(strrep(".", p5_start-1), p5bind, strrep(".", targetlen - p5_end)) %>% str_replace_all("\\(", "\\)")) %>%
    filter(str_detect(piRNAmatch, "\\(")) %>%
    mutate(contact = paste0(".",piRNAmatch,".",  targetmatch,"."),
           seq = paste0("A",piRNAseq, "A",targetseq, "A"))
table_RNA <- vienna_struct %>% select(contact, seq) %>% distinct()

cal_struct <- function(dat) {
  res <- dat %>% filter(group!=0) %>%
    group_by(group) %>% nest_legacy() %>%
    mutate(total = map(data,function(x) {x$num}),
           total2 = map(data,function(x) {x$num %>% diff()}),
           total3 = map(total2,function(x) {which(x>1)})) %>%
    select(-data) %>%
    filter(total3!="integer(0)")
  if (dim(res)[1]==0) {
    return(tibble(group = 0,before = 0,after = 0, anno = "perfect"))
  }else{
    res2 <- res %>% mutate(before = map2(total,total3, function(x,y) {x[1:y] %>% length()}),
                           after = map2(total, before, function(x,y) {length(x)-y})) %>%
      select(-contains("total")) %>%
      unnest_legacy() %>%
      mutate(anno = if_else((before==2 & after == 2) | (before==3 & after == 3), "mismatch",
                            if_else(before == 1 & after > 1, "target_bulge",
                                    if_else(before > 1 & after == 1, "piRNA_bulge",
                                            if_else(before>1&after>1,"interior","else")))))
    return(res2)
  }
}

table_RNA2 <- table_RNA %>%
  rowwise() %>%
  mutate(ct = makeCt(struct = contact, seq = seq) %>% ct2coord_up() %>% cal_struct() %>% .$anno %>% paste0(collapse = "@"))

tableRNA3 <- vienna_struct %>% left_join(table_RNA2, by = c("contact","seq")) %>%
  select(-contact, -seq)
tableRNA3 %>% write_tsv(paste0(outname,"_struct.tsv"))
