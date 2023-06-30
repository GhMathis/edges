library(tidyverse)
library(ggforce)
library(RColorBrewer)
library(cowplot)
library(purrr)
library(viridis)
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

viral_host_sharing_extended_df = read.csv("output/viral_host_sharing_extended_df.csv", header = T)


host_recap = viral_host_sharing_extended_df%>%
  select(c("HostOrder_zeta","HostOrder","host", "score_host_group","deepness","X"))%>%
  group_by(HostOrder,HostOrder_zeta, deepness)%>%
  summarise(mean_score = mean(score_host_group, na.rm=TRUE),
            sd_score = sd(score_host_group, na.rm=TRUE))

host_recap$deepness = as.factor(host_recap$deepness)

normalized_2=  function(x){
  return(x/sqrt(x**2+1000))
}
host_recap$mean_nrmlz = normalized_2(host_recap$mean_score)
host_recap %>% 
  filter(deepness == "Inf") %>% 
  ggplot()+
  geom_point(aes(HostOrder_zeta,HostOrder , col = normalized_2(mean_score)),size = 20,shape = 15)+
  scale_color_viridis()+
  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))

ggplot(host_recap)+
  geom_point(aes(HostOrder_zeta,deepness , col = mean_score),size = 10,shape = 15)+
  facet_wrap(~HostOrder)+

  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))


host_recap %>% 
  group_by(HostOrder) %>% 
  group_map(~ggplot(.)+geom_point( aes(HostOrder_zeta, deepness, color = mean_score),size = 10,shape = 15)+
              facet_wrap(~HostOrder)+
              scale_color_viridis(),.keep = T )%>%
  plot_grid(plotlist = ., align = 'hv', ncol = 3)

host_recap %>% 
  group_by(HostOrder) %>% 
  filter(deepness == "Inf") %>% 
  group_map(~ggplot(.)+
              geom_point( aes(HostOrder_zeta, HostOrder, group =HostOrder,
                                       col =mean_score),size = 10,shape = 15)+
              #facet_wrap(~HostOrder)+
              scale_color_viridis(),.keep = T)%>%
  plot_grid(plotlist = ., align = 'v', ncol = 1)

temp = viral_host_sharing_extended_df%>%
  select(c("HostOrder_zeta","HostOrder","host", "score_host_group","deepness","X"))%>%
  filter(HostOrder_zeta ==HostOrder )
  group_by(HostOrder,HostOrder_zeta, deepness)%>%
  summarise(mean_score = mean(score_host_group, na.rm=TRUE),
            sd_score = sd(score_host_group, na.rm=TRUE))
  
  
test = viral_host_sharing_extended_df%>%
    filter(HostOrder==HostOrder_zeta & deepness == "Inf")
ggplot(test)+
  geom_sina(aes(HostOrder_zeta, score_host_group))
  