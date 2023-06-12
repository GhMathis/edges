library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(car)
library(ggforce)
library(ggrepel)
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
host_subset2[host_subset2$X %in%c(640:1279),]$X == host_subset$X
host_subset2[host_subset2$X %in%c(640:1279),]$g_zeta_intra_order_h == host_subset$g_zeta_intra_order_h
temp = host_subset2[host_subset2$X %in%c(640:1279),]
host_subset = host_subset2
temp = host_subset
host_subset = read.csv("output/intra_vs_extra_host_subset.csv", header = T,
                                   stringsAsFactors = T)

host_subset$score_location = as.factor(host_subset$score_location)
levels(host_subset$score_location) = c("Extra Ordre", "Intra Ordre")
str(host_subset)
pal = brewer.pal(n =9, name = "OrRd")

as.factor(host_subset$HostOrder_zeta)
main_order = host_subset%>%
  count(HostOrder_zeta)%>%
  arrange(desc(n))%>%
  slice(1:10)%>%
  pull(HostOrder_zeta)
recap2 = host_subset%>%
  filter(HostOrder_zeta %in% unique(main_order))%>%
  group_by(score_location,HostOrder_zeta)%>%
  summarise(N    =length(score),
            mean_score = mean(score, na.rm=TRUE),
            sd_score = sd(score, na.rm=TRUE),
            se_score = sd_score / sqrt(N))
host_subset$HostOrder_zeta = ordered(host_subset$HostOrder_zeta, levels =unique(main_order)) 
host_subset%>%
  filter(HostOrder_zeta %in% unique(main_order))%>%
  #slice(1:2000)%>%
ggplot()+
  geom_hline(yintercept = 0, linetype =2)+
  geom_sina(aes(HostOrder_zeta, score, col = score_location),scale = "width", position=position_dodge(0.5),
            maxwidth = 0.4, cex =2, alpha = 0.5)+
  geom_errorbar(data = recap2, aes(x= HostOrder_zeta, group = score_location,
                                  ymin=mean_score-sd_score, ymax=mean_score+sd_score),
                cex =0.8, width=.3, position = position_dodge(0.5)) +
  geom_point(data = recap2, aes(x= HostOrder_zeta, y=mean_score,
                               group = score_location), cex =3,
             position=position_dodge(0.5)) +
  color_palette(c("#548C82",pal[8]))+
  #color_palette(pal[c(8,6,4)])+
  #scale_color_brewer("set1")+
  labs(col ="Calcule du partage viral", x = "Ordre", y = "z-score du partage viral potentiel")+
  main_theme+
  guides(colour = guide_legend(override.aes = list(alpha = 1, cex =3)))+
  theme(axis.text.x = element_text(angle =45, hjust = 1),
        legend.position = c(0.2, 0.85),
        legend.title = element_text(colour = "black", size=28),
        axis.title=element_text(size=32),
        legend.background = element_rect(fill = "transparent"),
        legend.box = "horizontal")


mod_v_sharing = glm(score~HostOrder_zeta*score_location, data = host_subset%>%
                     filter(HostOrder_zeta %in% unique(main_order)))

par(mfrow = c(2,2))
plot(mod_v_sharing)
hist(mod_v_sharing$residuals, breaks = 100)
summary(mod_v_sharing)
Anova(mod_v_sharing, type = 2)
TukeyHSD(aov(mod_v_sharing))

tablepresentation = recap2%>%
  select(-se_score)%>%
  mutate(mean_score = round(mean_score, 3),sd_score = round(sd_score, 3) )

  
###### Trefle vs Clover #######

## importance and communicability score for observed associassion (trefle and clover)
importance_df = read.csv("output/importance_df.csv")
summary(importance_df)
length(importance_df$status)
G_nrmlz_subset = importance_df %>%
  select(c(virus, host, X, G_pq_clover_nrmlz, G_pq_trefle_nrmlz))%>%
  pivot_longer(-c(virus, host, X),values_to ="G_nrmlz")
I_nrmlz_subset = importance_df %>%
  select(c(virus, host, X, importance_trefle_nrmlz, importance_clover_nrmlz))%>%
  pivot_longer(-c(virus, host, X),values_to ="I_nrmlz")

#### data whitout sapiens
## importance and communicability score for all associassion (trefle)
importance_trefle_nosapiens_df = read.csv("output/importance_trefle_nosapiens_df.csv")
importance_nosapiens_df = read.csv("output/importance_nosapiens_df.csv")


hist_G = ggplot(importance_df)+
  geom_density(aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz))+
  #geom_smooth(aes(G_pq_clover_nrmlz,G_pq_trefle_nrmlz), method = "lm", color = "black", linetype = 2, cex = 2)+
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  labs(fill ="Compte du nombre \n associations", x = "Communicabilité CLOVER", y= "Communicabilité TREFLE")+
  main_theme+
  theme(legend.title.align = 0.5,
        legend.direction = "vertical",
        legend.box.just = "center")+
  coord_fixed()


hist_G = ggplot(G_nrmlz_subset)+
  geom_density(aes(G_nrmlz, fill = name),alpha = 0.3, stat = "density")+
  fill_palette("pastel1")+
  main_theme+
  scale_fill_discrete(name = "Réseaux", labels = c("Clover", "Trefle"))+
  labs(x ="Communicabilité", y= "Densité")+
  theme(legend.position ="none")

hist_I = ggplot(I_nrmlz_subset)+
  geom_density(aes(I_nrmlz, fill = name),alpha = 0.3, stat = "density")+
  fill_palette("pastel1")+
  scale_fill_discrete(name = "Réseaux", labels = c("Clover", "Trefle"))+
  labs(x ="Importance", y= "Densité") +
  main_theme+
  theme(legend.position = c(0.7, 0.4),
        legend.text = element_text(colour = "black", size=12))
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
subset_gg_virusmax[which(subset_gg_virusmax$host == "Homo sapiens"),]
str(subset_gg_rest)
subset_gg_rest = importance_df%>%
  filter(!(host %in% host_max[1:10]),!(virus %in% virus_max[1:2]))



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
  
subset_gg_I2 = subset(importance_df, importance_trefle_nrmlz > 0.01,
                        importance_clover_nrmlz > 0.0405)

str(subset_gg_rest)
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
 

recap_trefle_clover = importance_df %>%
  summarise(across(c(importance_trefle_nrmlz, importance_clover_nrmlz,
                     G_pq_trefle_nrmlz, G_pq_clover_nrmlz),
                   list(mean = mean, sd = sd)))%>%
  pivot_longer(everything())
# SE
(importance_df %>%
  summarise(across(c(importance_trefle_nrmlz, importance_clover_nrmlz,
                     G_pq_trefle_nrmlz, G_pq_clover_nrmlz),
                   list(mean = mean, sd = sd)))%>%
  select(ends_with("_sd")))/sqrt(nrow(importance_df))
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
#virus and host centrality
importance_df = importance_df%>%
  full_join(importance_df%>%
              count(virus),by = "virus" )%>%
  full_join(importance_df%>%
              count(host),by = "host" )
names(importance_df)[20:21] = c("c_virus", "c_host")
gam_I_N = mgcv::gam(importance_trefle_nrmlz ~ s(importance_clover_nrmlz,bs = "cs")+c_host*c_virus,
                  data = importance_df)
summary(gam_I_N)
plot(gam_I_N)
plot(gam_I)
AIC(gam_I)
AIC(gam_I_N)
summary(gam_I_N)
summary(gam_I)

plot(gam_I$fitted.values,gam_I$residuals)
# Add a line at y = 0
abline(h = 0, col = "red", lty = 2)
hist(gam_I$residuals~gam_I$fitted.values)
preds = data.frame(importance_clover_nrmlz=seq(from=min(importance_df$importance_clover_nrmlz), 
    to=max(importance_df$importance_clover_nrmlz),length.out = 1000),
    c_host = sample(min(importance_df$c_host):max(importance_df$c_host),1000, replace = T ),
    c_virus= sample(min(importance_df$c_virus):max(importance_df$c_virus),1000, replace = T  ))

preds = cbind(preds,as_tibble(
  predict(gam_I_N, newdata = preds, se.fit=TRUE)))


str(preds)
gg_trefle = ggplot(importance_df)+
  geom_point(aes(G_pq_trefle_nrmlz, importance_trefle_nrmlz))+
  geom_abline(intercept= 0, slope = 1, col = "blue")+
  labs(x ="G trefle", y= "I trefle")+
  main_theme

### Z value distribution comparaison

sum(subset_gg_virusmax$G_pq_trefle_nrmlz<subset_gg_virusmax$G_pq_clover_nrmlz)
sum(importance_df$G_pq_trefle_nrmlz<importance_df$G_pq_clover_nrmlz)
sum(importance_df$G_pq_trefle_nrmlz<importance_df$G_pq_clover_nrmlz
  & importance_df$host == "Homo sapiens")
sum(importance_df$G_pq_trefle_nrmlz>=importance_df$G_pq_clover_nrmlz
    & importance_df$virus == "West Nile virus")
sum(importance_df$G_pq_trefle_nrmlz>importance_df$G_pq_clover_nrmlz
    & importance_df$virus == "Rabies lyssavirus")

sum(importance_df$host == "Homo sapiens")
sum( importance_df$virus == "West Nile virus")
1-NA
mean(na.rm = T)
nrow(importance_df)
( mean(importance_df$G_pq_trefle_nrmlz, na.rm = )- 
    mean(importance_df$G_pq_clover_nrmlz))/sqrt(
      (sd(importance_df$G_pq_clover_nrmlz)/sqrt(nrow(importance_df)))^2 +
       (sd(importance_df$G_pq_trefle_nrmlz)/sqrt(nrow(importance_df)))^2
    )

(mean(importance_df$importance_trefle_nrmlz)- 
    mean(importance_df$importance_clover_nrmlz))/sqrt(
      (sd(importance_df$importance_trefle_nrmlz)/sqrt(nrow(importance_df)))^2 +
       (sd(importance_df$importance_clover_nrmlz)/sqrt(nrow(importance_df)))^2
    )


###### Importance vs chance of association (trefle only) #######
imputed_association = read.csv("data/imputed_associations.csv")
#### data with sapiens
## importance and communicability score for imputed associassion (trefle only)
importance_trefle_df = read.csv("output/importance_unshared_df.csv")
importance_trefle_df$status = "imputed associations"

## importance and communicability score for observed associassion (trefle and clover)
importance_df = read.csv("output/importance_df.csv")
importance_df$status = "observed associations"

#### data whitout sapiens
# data whitout sapiens

importance_trefle_nosapiens_df = read.csv("output/importance_trefle_nosapiens_df.csv")
importance_nosapiens_df = read.csv("output/importance_nosapiens_df.csv")

v = importance_trefle_nosapiens_df$virus
h = importance_trefle_nosapiens_df$host

ID = which(paste(v,h) %in% paste(imputed_association$virus, imputed_association$host))
importance_trefle_nosapiens_df$status[ID] = "imputed associations"
importance_trefle_nosapiens_df$status[-ID] = "observed associations"
importance_nosapiens_df[ID,]


# joint to have association proba and importance in same the df
str(imputed_association)
imputed_association = imputed_association%>%
  full_join(importance_trefle_df,by =c("host", "virus"))
str(imputed_association)

# Imputed + observed associations
importance_trefle_df = rbind(importance_trefle_df, importance_df[,names(importance_df) %in% names(importance_trefle_df)])

str(importance_trefle_df)

ggplot(imputed_association)+
  geom_bin_2d(aes(G_pq_trefle_nrmlz,evidence), bins = 100)+
  geom_smooth(aes(G_pq_trefle_nrmlz,evidence), method = "lm", color = "black", linetype = 2, cex = 2)+
  scale_y_log10()+
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  labs(fill ="Compte du nombre \n associations", x = "Communicabilité", y= "evidence (échelle log)")+
  main_theme+
  theme(legend.title.align = 0.5,
        legend.direction = "vertical",
        legend.box.just = "center")

######
# Imputed + observed associations

str(importance_trefle_df)
# imputed + observed and ordering to have coresspndig col names
importance_trefle_df = rbind(importance_trefle_df,
                             importance_df[names(importance_df) %in% names(importance_trefle_df)])
gg_sap = ggplot(importance_trefle_df)+
  #geom_bin_2d(aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz))+
  geom_point(aes(G_pq_trefle_nrmlz, importance_trefle_nrmlz, col= status), alpha = 0.1,cex =2)+
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  # geom_text_repel(data=subset_gg_G2,
  #                 aes(G_pq_clover_nrmlz,G_pq_trefle_nrmlz, label = paste(host, "(", virus, ")", sep = "")),
  #                 size=4,
  #                 box.padding = unit(0.5, "lines"),
  #                 max.overlaps =20)+
  coord_cartesian(xlim = c(0,0.77))+
  labs(x ="G trefle", y= "I trefle")+
  main_theme


gg_sap_nosap = ggplot(importance_trefle_nosapiens_df)+
  #geom_bin_2d(aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz))+
  geom_point(aes(G_pq_trefle_nrmlz, importance_trefle_nrmlz, col= status), alpha = 0.1,cex =2)+
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  # geom_text_repel(data=subset_gg_G2,
  #                 aes(G_pq_clover_nrmlz,G_pq_trefle_nrmlz, label = paste(host, "(", virus, ")", sep = "")),
  #                 size=4,
  #                 box.padding = unit(0.5, "lines"),
  #                 max.overlaps =20)+
  coord_cartesian(xlim = c(0,0.77))+
  labs(x ="G trefle", y= "I trefle")+
  main_theme

#G
importance_trefle_nosapiens_df%>%
  arrange(desc(G_pq_trefle_nrmlz))%>%
  filter(status == "imputed associations")%>%
  slice(1:20)%>%
  select(host,virus,status)
importance_trefle_df%>%
  arrange(desc(G_pq_trefle_nrmlz))%>%
  filter(status == "imputed associations")%>%
  slice(1:20)%>%
  select(host,virus,status)
#I
importance_trefle_df%>%
  arrange(desc(importance_trefle_nrmlz))%>%
  #filter(status == "imputed associations")%>%
  slice(1:20)%>%
  select(host,virus,status)
importance_trefle_nosapiens_df%>%
  arrange(desc(importance_trefle_nrmlz))%>%
  #filter(status == "imputed associations")%>%
  slice(1:20)%>%
  select(host,virus,status)


load("output/G_trefle.Rdata")
clover = read.csv("data/clover.csv", header = T,
                  stringsAsFactors = T)
ID_host= which(colnames(G_trefle) %in% unique(clover$Host))
max(G_trefle[ID_host,ID_host])


# Extract the maximum values and their corresponding row and column names
max_values <- mat[max_indices]
row_names <- row.names(mat)[max_indices[, 1]]
col_names <- colnames(mat)[max_indices[, 2]]


###### Network

source("R/functions.R")
clover = read.csv("data/clover.csv")
trefle = read.csv("data/trefle.csv")
uni_ntw_clover = matrix.associations.uni(Virus = clover$Virus, Host = clover$Host)
uni_ntw_trefle = matrix.associations.uni(Virus = trefle$virus, Host = trefle$host)

grh_clover = graph_from_adjacency_matrix(uni_ntw_clover,mode = "undirected")
grh_trefle = graph_from_adjacency_matrix(uni_ntw_trefle,mode = "undirected")

str(summary(importance_df$importance_trefle_nrmlz))

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


edges_t <- degree(grh_trefle)
colorvector = c()
for(i in 1:length(edges_t)){
  alph = edges_t[i]
    colorvector = c(colorvector, adjustcolor("gray80",alpha =alph))
}
edges.color = colorvector



I_nomr_clvr = importance_df$importance_clover_nrmlz
edge_color_share2 = ifelse(I_nomr_clvr>=max(I_nomr_trfl) - max(I_nomr_trfl)*0.01 , "blue",
                           ifelse(I_nomr_clvr==0, brewer.pal(n=9,name ="YlOrRd" )[2],
                                  ifelse(I_nomr_clvr<max(I_nomr_trfl) - max(I_nomr_trfl)*0.2, brewer.pal(n=9,name ="YlOrRd" )[6],
                                         ifelse(I_nomr_clvr<max(I_nomr_trfl) - max(I_nomr_trfl)*0.1, brewer.pal(n=9,name ="YlOrRd" )[7],
                                                brewer.pal(n=9,name ="YlOrRd" )[9])
                                  )))

degree_c <- degree(grh_clover)
plot(grh_clover, vertex.size=log(degree_c),vertex.label=NA,
     vertex.label.cex = 2,edge.width = 1,
     edge.label=NA,edge.label.cex = 2,
     edge.color = edge_color_share, layout = coords)

##### Edges value from clover

plot(grh_clover, vertex.size=1,vertex.label=NA,
     vertex.label.cex = 2,edge.width = 8,
     edge.label=NA,edge.label.cex = 2,
     edge.color = edge_color_share2, layout = lay)