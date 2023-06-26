library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(car)
library(ggforce)
library(ggrepel)
library(igraph)
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

###### Viral sharing #######
viral_host_sharing_df = read.csv("output/viral_host_sharing_dG_df.csv", header = T,
                                   stringsAsFactors = T)
host_subset = viral_host_sharing_df%>%
  select(c(HostOrder_zeta, host, score_intra_h, score_extra_h))%>%
  pivot_longer(c(score_intra_h,score_extra_h), names_to ="score_location", values_to = "score")

host_subset$score_location = as.factor(host_subset$score_location)
levels(host_subset$score_location) = c("Extra Ordre", "Intra Ordre")
str(host_subset)

### arrange orders by numbers of associations
as.factor(host_subset$HostOrder_zeta)
main_order = host_subset%>%
  count(HostOrder_zeta)%>%
  arrange(desc(n))%>%
  slice(1:10)%>%
  pull(HostOrder_zeta)
recap = host_subset%>%
  filter(HostOrder_zeta %in% unique(main_order))%>%
  group_by(score_location,HostOrder_zeta)%>%
  summarise(N    =length(score),
            mean_score = mean(score, na.rm=TRUE),
            sd_score = sd(score, na.rm=TRUE),
            se_score = sd_score / sqrt(N))
host_subset$HostOrder_zeta = ordered(host_subset$HostOrder_zeta, levels =unique(main_order)) 

host_subset%>%
  filter(HostOrder_zeta %in% unique(main_order))%>%
ggplot()+
  geom_hline(yintercept = 0, linetype =2)+
  geom_sina(aes(HostOrder_zeta, score, col = score_location),scale = "width", position=position_dodge(0.5),
            maxwidth = 0.4, cex =2, alpha = 0.5)+
  geom_errorbar(data = recap, aes(x= HostOrder_zeta, group = score_location,
                                  ymin=mean_score-sd_score, ymax=mean_score+sd_score),
                cex =0.8, width=.3, position = position_dodge(0.5)) +
  geom_point(data = recap, aes(x= HostOrder_zeta, y=mean_score,
                               group = score_location), cex =3,
             position=position_dodge(0.5)) +
  color_palette(c("#548C82",pal[8]))+
  labs(col ="Calcule du partage viral", x = "Ordre", y = "z-score du partage viral potentiel")+
  main_theme+
  guides(colour = guide_legend(override.aes = list(alpha = 1, cex =3)))+
  theme(axis.text.x = element_text(angle =45, hjust = 1),
        legend.position = c(0.2, 0.85),
        legend.title = element_text(colour = "black", size=30),
        axis.title=element_text(size=32),
        axis.text = element_text(colour = "black", size=26),
        legend.text = element_text(colour = "black", size=26),
        legend.background = element_rect(fill = "transparent"),
        legend.box = "horizontal")

### glm
mod_v_sharing = glm(score~HostOrder_zeta*score_location, data = host_subset%>%
                     filter(HostOrder_zeta %in% unique(main_order)))
par(mfrow = c(2,2))
plot(mod_v_sharing)
hist(mod_v_sharing$residuals, breaks = 100)
summary(mod_v_sharing)
Anova(mod_v_sharing, type = 2)


### table of N_asso, means and sd
tablepresentation = recap%>%
  select(-se_score)%>%
  mutate(mean_score = round(mean_score, 3),sd_score = round(sd_score, 3) )

  
###### Trefle vs Clover #######

## importance and communicability score for observed associassion (trefle and clover)
importance_df = read.csv("output/importance_df.csv")

str(importance_df)
G_nrmlz_subset = importance_df %>%
  select(c(virus, host, X, G_pq_clover_nrmlz, G_pq_trefle_nrmlz,))%>%
  pivot_longer(-c(virus, host, X),values_to ="G_nrmlz")
I_nrmlz_subset = importance_df %>%
  select(c(virus, host, X, importance_trefle_nrmlz, importance_clover_nrmlz))%>%
  pivot_longer(-c(virus, host, X),values_to ="I_nrmlz")

### density communicability
hist_G = ggplot(G_nrmlz_subset)+
  geom_density(aes(G_nrmlz, fill = name),alpha = 0.3, stat = "density")+
  fill_palette("pastel1")+
  main_theme+
  scale_fill_discrete(name = "Réseaux", labels = c("Clover", "Trefle"))+
  labs(x ="Communicabilité", y= "Densité")+
  theme(legend.position ="none")

### density importance
hist_I = ggplot(I_nrmlz_subset)+
  geom_density(aes(I_nrmlz, fill = name),alpha = 0.3, stat = "density")+
  fill_palette("pastel1")+
  scale_fill_discrete(name = "Réseaux", labels = c("Clover", "Trefle"))+
  labs(x ="Importance", y= "Densité") +
  main_theme+
  theme(legend.position = c(0.7, 0.4),
        legend.text = element_text(colour = "black", size=12))

### Virus and host with most association (for the legend)
importance_df$host = as.factor(importance_df$host)
host_max = importance_df%>%
  count(host)%>%
  arrange(desc(n))%>%pull(host)
importance_df$host <- ordered(importance_df$host , levels =unique(host_max))
subset_gg_hostmax = importance_df%>%
  filter(host %in% host_max[1:10])
virus_max = importance_df%>%
  count(virus)%>%
  arrange(desc(n))%>%
  #slice(1:2)
  pull(virus)
subset_gg_virusmax = importance_df%>%
  filter(virus %in% virus_max[1:2])

#### Virus and host with less associations
subset_gg_rest = importance_df%>%
  filter(!(host %in% host_max[1:10]),!(virus %in% virus_max[1:2]))


##### scatter plot communicability
subset_gg_G2 = subset(importance_df, G_pq_trefle_nrmlz > 0.86|
                        G_pq_clover_nrmlz > 0.3)
gg_G2  = ggplot(subset_gg_rest)+
  #geom_bin_2d(aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz))+
  geom_abline(intercept= 0, slope = 1, col = "red", linetype =2, cex =1)+
  geom_point(aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz),
             alpha = 0.1,cex =3)+
  geom_point(data=subset_gg_hostmax,
             aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz,col = host), 
             alpha = 0.4, cex =3)+
  geom_point(data=subset_gg_virusmax,
             aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz,shape = virus),
             alpha =0.6, cex =3)+
  scale_shape_manual(values = c(5,6))+
  
  #scale_color_manual(values = brewer.pal(10, "Set3")) +
  
  geom_text_repel(data=subset_gg_G2,
                  aes(G_pq_clover_nrmlz,G_pq_trefle_nrmlz, label = paste(host, "(", virus, ")", sep = "")),
                  size=4,
                  box.padding = unit(1, "lines"),
                  max.overlaps =20)+
  labs(x ="Communicabilité Clover", y= "Communicabilité Trefle",
       col = "Hôtes", shape = "Virus")+
  main_theme+
  theme(legend.position ="none")


##### scatter plot importance  
subset_gg_I2 = subset(importance_df, importance_trefle_nrmlz > 0.01,
                        importance_clover_nrmlz > 0.0405)
gg_I2 = ggplot(subset_gg_rest)+
  geom_abline(intercept= 0, slope = 1, col = "red", linetype =2, cex =1)+
  geom_point(aes(importance_clover_nrmlz, importance_trefle_nrmlz),
             alpha = 0.1,cex =3)+
  geom_point(data=subset_gg_hostmax,
             aes(importance_clover_nrmlz, importance_trefle_nrmlz,col = host),
             alpha = 0.4, cex =3)+
  geom_point(data=subset_gg_virusmax,
             aes(importance_clover_nrmlz, importance_trefle_nrmlz,shape = virus),
             alpha =0.6, cex =3)+
  scale_shape_manual(values = c(5,6))+
  geom_text_repel(data=subset_gg_I2,
                  aes(importance_clover_nrmlz,importance_trefle_nrmlz, label = paste(host, "(", virus, ")", sep = "")),
                  size=4,
                  box.padding = unit(0.5, "lines"),
                  max.overlaps = 20)+

  labs(x ="Importance Clover", y= "Importance Trefle",
       col = "Hôtes", shape = "Virus")+
  #geom_line(data =preds, aes(x=importance_clover_nrmlz, y= fit))+ 
  geom_smooth(data = importance_df, aes(importance_clover_nrmlz, importance_trefle_nrmlz),
              cex =2)+
  main_theme+
  guides(colour = guide_legend(override.aes = list(alpha = 0.8, cex =3),ncol = 2))+
  theme(legend.position = c(0.62, 0.2),
        legend.text = element_text(colour = "black", size=12),
        legend.background = element_rect(fill = "transparent"),
        legend.box = "horizontal")
  
  

ggarrange(hist_G,hist_I,gg_G2,gg_I2,labels =c("A","B","C","D"),
          label.x = 0.2, ncol =2,nrow = 2, align = "v",
          font.label = list(size = 20, face = "bold"))
 
### mean and sd
recap_trefle_clover = importance_df %>%
  summarise(across(c(importance_trefle_nrmlz, importance_clover_nrmlz,
                     G_pq_trefle_nrmlz, G_pq_clover_nrmlz),
                   list(mean = mean, sd = sd)))%>%
  pivot_longer(everything())
# model Itrefle~log(clover)

str(importance_df)

help(summary.gam)
mod_I = lm(importance_trefle_nrmlz ~ importance_clover_nrmlz,
                  data = importance_df)
mod2_I = mgcv::gam(importance_trefle_nrmlz ~ poly(importance_clover_nrmlz,2),
           data = importance_df)
gam_I2 = mgcv::gam(importance_trefle_nrmlz ~ s(importance_clover_nrmlz,bs = "cr"),
                  data = importance_df)
gam_I = mgcv::gam(importance_trefle_nrmlz ~ s(importance_clover_nrmlz,bs = "cs"),
            data = importance_df)
AIC(gam_I2)
AIC(gam_I)
plot(gam_I2)
plot(gam_I)

###### Importance vs chance of association (trefle only) #######
imputed_association = read.csv("data/imputed_associations.csv")
#### data with sapiens
## importance and communicability score for imputed associassion (trefle only)
importance_trefle_df = read.csv("output/importance_unshared_df.csv")
importance_trefle_df$status = "imputed associations"

## importance and communicability score for observed associassion (trefle and clover)
importance_df = read.csv("output/importance_df.csv")
importance_df$status = "observed associations"


###### Network of imputed association

source("R/functions.R")
clover = read.csv("data/clover.csv")
trefle = read.csv("data/trefle.csv")
uni_ntw_clover = matrix.associations.uni(Virus = clover$Virus, Host = clover$Host)
uni_ntw_trefle = matrix.associations.uni(Virus = trefle$virus, Host = trefle$host)

grh_clover = graph_from_adjacency_matrix(uni_ntw_clover,mode = "undirected")
grh_trefle = graph_from_adjacency_matrix(uni_ntw_trefle,mode = "undirected")



ID_share_edges = which(paste(as_edgelist(grh_trefle, names=T)[,1],
                             as_edgelist(grh_trefle, names=T)[,2]) %in% 
                         paste(as_edgelist(grh_clover, names=T)[,1],
                               as_edgelist(grh_clover, names=T)[,2]))

I_nomr_trfl = importance_df$importance_trefle_nrmlz
edge_color_share = ifelse(I_nomr_trfl>=max(I_nomr_trfl) - max(I_nomr_trfl)*0.01 , "blue",
                          ifelse(importance_df$importance_trefle_nrmlz==0, brewer.pal(n=9,name ="YlOrRd" )[2],
                                 ifelse(I_nomr_trfl<max(I_nomr_trfl) - max(I_nomr_trfl)*0.2, brewer.pal(n=9,name ="YlOrRd" )[6],
                                        ifelse(I_nomr_trfl<max(I_nomr_trfl) - max(I_nomr_trfl)*0.1, brewer.pal(n=9,name ="YlOrRd" )[7],
                                               brewer.pal(n=9,name ="YlOrRd" )[9])
                                 )))
edge_color_trefle = rep("gray80",length(E(grh_trefle)))
edge_color_trefle[ID_share_edges] <- edge_color_share

edge_width_trefle = rep(8,length(E(grh_trefle)))
edge_width_trefle[-ID_share_edges] =NA

##### Edges value from trefle
set.seed(1)
coords = layout_(grh_clover,with_fr(niter = 10000),x = component_wise())
coordstr = layout_(grh_trefle,with_fr(niter = 10000),x = component_wise())
degree_t <- degree(grh_trefle)
colorvector = c()
for(i in 1:length(degree_t)){
  alph = c(degree_t/max(degree_t))[i]
  if(attr(alph,"name" ) %in% c(unique(trefle$virus),unique(clover$Virus))){
    colorvector = c(colorvector, adjustcolor("#FFDFA0",alpha =alph))
  }else{
    colorvector = c(colorvector, adjustcolor("#E3C6FF",alpha =alph))
  }
}

vertex.color = colorvector
I_nomr_trfl = importance_df$importance_trefle_nrmlz
edge_color_share = ifelse(I_nomr_trfl>=max(I_nomr_trfl) - max(I_nomr_trfl)*0.05 ,  adjustcolor("blue",1),
                          ifelse(importance_df$importance_trefle_nrmlz==0, adjustcolor("#FFEDA0", 0.1),
                                 ifelse(I_nomr_trfl<max(I_nomr_trfl) - max(I_nomr_trfl)*0.2, adjustcolor("#FC4E2A", 0.2),
                                        ifelse(I_nomr_trfl<max(I_nomr_trfl) - max(I_nomr_trfl)*0.1, adjustcolor("#E31A1C",0.3),
                                               adjustcolor("#800026", 0.8))
                                 )))

plot(grh_trefle, vertex.size=sqrt(log(degree_t)),vertex.label=NA,
     vertex.color = vertex.color,
     vertex.label.cex = 2,edge.width = 1,
     edge.label=NA,edge.label.cex = 2,
     edge.color = edge_color_share, layout = coordstr)
