library(colorRamps)
library(tidyverse)
library(ggforce)
source("R/functions.R")
main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20),
        legend.title.align=0.5,
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=28))

pred_ntw = readRDS("data/PredictedNetwork.rds")
host_data = read.csv("data/Hosts.csv", header = T)

str(pred_ntw)
str(host_data)
mean_sharing = mean(pred_ntw)
sd_sharing = sd(pred_ntw)
score_host_group = c()
HostOrder1 = c()
HostOrder2 = c()
for(HO in unique(host_data$hOrder)){
  for(HO2 in unique(host_data$hOrder)){
  #
  host_group = host_data$hHostNameFinal[host_data$hOrder == HO]
  host_group2 = host_data$hHostNameFinal[host_data$hOrder == HO2]
  
  ID_host_group =  which(colnames(pred_ntw) %in% host_group)
  ID_host_group2 =  which(rownames(pred_ntw) %in% host_group2)
  
  mean_host_group = mean(pred_ntw[ID_host_group2, ID_host_group])
  
  score_host_group = c(score_host_group, (mean_host_group- 
                                            mean_sharing)/sd_sharing)
  HostOrder1 =c(HostOrder1, HO)
  HostOrder2 = c(HostOrder2, HO2)
  }
}
df_sharing = data.frame(HostOrder1 = HostOrder1, HostOrder2 = HostOrder2,
                        score = score_host_group)
df_sharing = na.omit(df_sharing)
df_sharing$score_nrmlz = normalized_2(df_sharing$score)
df_sharing = df_sharing%>%
  group_by(HostOrder2)%>%
  mutate(score_nrmlz_grp = normalized_2(score))

ggplot(df_sharing)+
  geom_point(aes(HostOrder1,HostOrder2 , col = score_nrmlz),size = 10,shape = 15)+
  scale_colour_gradientn(colours = rev(green2red(500)),limits=c(-1, 1))+
  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(df_sharing)+
  geom_point(aes(HostOrder1,HostOrder2 , col = score_nrmlz_grp),size = 10,shape = 15)+
  scale_colour_gradientn(colours = rev(green2red(500)),limits=c(-1, 1))+
  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))
