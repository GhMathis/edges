library(tidyverse)
library(ggforce)
library(RColorBrewer)
library(cowplot)
library(purrr)
library(colorRamps)
library(lattice)
library(viridis)
library(ggforce)
library(reshape2)
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

str(viral_host_sharing_extended_df)
host_recap = viral_host_sharing_extended_df%>%
  filter(HostOrder_zeta %in% HostOrder)%>%
  select(c("HostOrder_zeta","HostOrder","host", "score_host_group","deepness","X"))%>%
  group_by(HostOrder,HostOrder_zeta, deepness)%>%
  summarise(mean_score = mean(score_host_group, na.rm=TRUE),
            sd_score = sd(score_host_group, na.rm=TRUE))

host_recap$deepness = as.factor(host_recap$deepness)

viral_host_sharing_extended_df%>%
  filter( deepness == "Inf")%>%
  ggplot()+
  geom_sina(aes(HostOrder_zeta, score_host_group, col = HostOrder))
str(host_recap)

host_recap = host_recap%>%
  group_by(HostOrder, deepness)%>%
  mutate(mean_nrmlz = normalized_2(mean_score))



host_recap_inf = host_recap %>% 
  filter(deepness == "Inf")
names(airquality) <- tolower(names(airquality))
aqm <- melt(airquality, id=c("month", "day"), na.rm=TRUE)

host_recap_mtx_nrmlz = acast(host_recap_inf%>%
                         ungroup()%>%
                        select(HostOrder, HostOrder_zeta, mean_nrmlz),
                       HostOrder~HostOrder_zeta)

color_scale1 = make.color.scale(host_recap_mtx_nrmlz)
levelplot(host_recap_mtx_nrmlz,
          scales = list( x=list(rot=90)),
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal,
          xlab = "spillover vulnerability", ylab = "spillover potetial" )

host_recap_mtx = acast(host_recap_inf%>%
                         ungroup()%>%
                         select(HostOrder, HostOrder_zeta, mean_score),
                       HostOrder~HostOrder_zeta)
color_scale2 = make.color.scale(host_recap_mtx)
levelplot(host_recap_mtx,
          scales = list( x=list(rot=90)),
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal,
          xlab = "spillover vulnerability", ylab = "spillover potetial" )


host_recap %>% 
  filter(deepness == "Inf") %>% 
  ggplot()+
  geom_point(aes(HostOrder_zeta,HostOrder , col = mean_nrmlz),size = 10,shape = 15)+
  scale_colour_gradientn(colours = rev(green2red(500)),limits=c(-1, 1) )+
  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))
host_recap %>% 
  filter(deepness == "Inf") %>% 
  ggplot()+
  geom_point(aes(HostOrder_zeta,HostOrder , col = mean_score),size = 10,shape = 15)+
  scale_colour_gradientn(colours = rev(green2red(500)), limits=c(-3, 3))+
  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(host_recap)+
  geom_point(aes(HostOrder_zeta,deepness , col = mean_nrmlz),size = 10,shape = 15)+
  facet_wrap(~HostOrder)+
  scale_colour_gradientn(colours = rev(green2red(500)),limits=c(-1, 1) )+
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
  