library(tidyverse)
library(lattice)
library(igraph)
library(corrplot)
library(colorRamps)
library(ggpubr)
library(RColorBrewer)
library(scico)
library(reshape2)
library(ggforce)

library(ape)
library(PVR)
library(phylogram)
source("R/functions.R")

main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=26),
        axis.text.y = element_text(colour = "black", size=26),
        legend.title = element_text(colour = "black", size=30),
        legend.title.align=0.5,
        legend.text = element_text(colour = "black", size=24),
        axis.title=element_text(size=28))

## aproximalty 0.28 To of storage needed
clover = read.csv("data/clover.csv", header = T)
trefle = read.csv("data/trefle.csv", header = T)
str(clover)
# trefle = trefle%>%
#   filter(host != "Homo sapiens")

##### compute adjacency matrix 
uni_ntw_trefle = matrix.associations.uni(Virus = trefle$virus, Host = trefle$host)
uni_ntw_clover = matrix.associations.uni(Virus = clover$Virus, Host = clover$Host)
##### compute communicability matrix
G_trefle = communicability(uni_ntw_trefle)
G_cluster = clustering.func(uni_ntw_trefle)
G_cluster = normalized_2(G_cluster)
G_cluster_host = G_cluster[rownames(G_cluster)%in% unique(trefle$host),
                           colnames(G_cluster)%in% unique(trefle$host)]
##### Load and compute phylo dist between host
mam_supertree <- read.tree("data/supertree_mammals.tree")
mam_supertree$tip.label = sub("_", " ", mam_supertree$tip.label)
ST_drop <- setdiff(mam_supertree$tip.label, unique(trefle$host)) #sp to drop in the tree (not presnete in trefle)
mam_supertree2 <- drop.tip(mam_supertree,
                           which(mam_supertree$tip.label %in% ST_drop)) #new tree
which(unique(trefle$host) == "Abrothrix longipilis")
##### Arrange host data of phylo dist and clustering  

phylodist <- cophenetic(mam_supertree2)
phylodist_df = melt(phylodist, value.name = "phylodist")

phylodist_df$ID_hosts = paste(phylodist_df$Var1, phylodist_df$Var2)
phylodist_df = phylodist_df[,!(colnames(phylodist_df) %in% c("Var1", "Var2"))]
str(phylodist_df)

df_clustering_host = melt(G_cluster_host[rownames(G_cluster_host) %in% mam_supertree$tip.label,
                                         colnames(G_cluster_host) %in% mam_supertree$tip.label],value.name = "clustering")

df_clustering_host$ID_hosts = paste(df_clustering_host$Var1, df_clustering_host$Var2)
df_clustering_host = df_clustering_host%>%
  full_join(phylodist_df, by = c("ID_hosts"))

#Communicability 
df_G_host = melt(G_trefle[rownames(G_trefle) %in% mam_supertree$tip.label,
                                         colnames(G_trefle) %in% mam_supertree$tip.label],value.name = "G")
df_G_host$ID_hosts = paste(df_G_host$Var1, df_G_host$Var2)
df_G_host = df_G_host[,!(colnames(df_G_host) %in% c("Var1", "Var2"))]
df_clustering_host = df_clustering_host%>%
  full_join(df_G_host, by = c("ID_hosts"))
# Add 10 most represented order to the df 
main_h_order = clover%>%
  count(HostOrder)%>%
  arrange(desc(n))%>%
  slice(1:10)
str(df_clustering_host)
df_clustering_host$HostOrder1 = "_"
for(h in unique(df_clustering_host$Var1)){
    HO = clover$HostOrder[which(clover$Host == h)[1]]
    host_indexs = which(df_clustering_host$Var1 == h)
    if (HO %in% main_h_order$HostOrder){
      df_clustering_host$HostOrder1[host_indexs] = as.character(HO)
    }else{
      df_clustering_host$HostOrder1[host_indexs] = "Other"
    }
}
df_clustering_host$HostOrder2 = "_"
for(h in unique(df_clustering_host$Var2)){
  HO = clover$HostOrder[which(clover$Host == h)[1]]
  host_indexs = which(df_clustering_host$Var2 == h)
  if (HO %in% main_h_order$HostOrder){
    df_clustering_host$HostOrder2[host_indexs] = as.character(HO)
  }else{
    df_clustering_host$HostOrder2[host_indexs] = "Other"
  }
}

str(df_clustering_host)
ggplot(df_clustering_host)+
  geom_bin_2d(aes(phylodist, clustering),bins = 50)+
  geom_smooth(aes(phylodist, clustering, col = HostOrder1), method ="gam")+
  color_palette("set3")+
  scale_x_log10()+
  facet_wrap(~HostOrder1)

# ggplot(df_clustering_host)+
#   geom_point(aes(clustering, phylodist, col = HostOrder2), alpha = 0.4)+
#   geom_smooth(aes(clustering, phylodist, col = HostOrder1))+
#   facet_wrap(~HostOrder1) 

df_clustering_host%>%
  filter(Var1 == "Homo sapiens")%>%
ggplot()+
  geom_point(aes(phylodist, clustering, col = HostOrder2), size = 5, alpha = 0.8 )+
  geom_smooth(aes(phylodist, clustering), method = "lm", cex = 1.5, linetype = 2)+
  geom_hline(yintercept  = 0, col = "black", cex = 1.25, linetype =2)+
  scale_x_log10()+
  scale_colour_brewer(palette = "Set3")+
  labs(x = "Phylogenetic distance (log scale)", y = "viral similarity", col = "Host Order of the 2nd part")+
  main_theme
str(df_clustering_host)
unique(df_clustering_host$HostOrder2)
str(recap_primate)
recap_primate = df_clustering_host%>%
  filter(HostOrder1 == "Primates")%>%
  group_by(HostOrder2) %>%
  summarise(N = length(clustering),
            mean_clustering = mean(clustering, na.rm=TRUE),
            sd_clustering = sd(clustering, na.rm=TRUE),
            se_clustering = clustering / sqrt(N),.groups = "keep")
df_clustering_host%>%
  filter(HostOrder1 == "Primates")%>%
ggplot()+
  geom_sina(aes(HostOrder2, clustering))+
  geom_hline(yintercept  = 0, col = "red", cex = 1.5, linetype =2)+
  scale_colour_brewer(palette = "Set3")+
  labs(x = "Host Order 2", y = "viral similarity")+
  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))


g1 = df_clustering_host%>%
  filter(Var1 == "Homo sapiens")%>%
  ggplot()+
  geom_boxplot(aes(HostOrder2, G))+
  scale_colour_brewer(palette = "Set3")+
  labs(x = "Host Order 2", y = "Communicability")+
  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))
g2 = df_clustering_host%>%
  filter(Var1 == "Homo sapiens")%>%
  ggplot()+
  geom_boxplot(aes(HostOrder2, clustering))+
  scale_colour_brewer(palette = "Set3")+
  labs(x = "Host Order 2", y = "viral spectrum similarity (clustering)")+
  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))

ggarrange(g1, g2)

ggplot(df_clustering_host)+
  geom_boxplot(aes(HostOrder2, G))+
  scale_colour_brewer(palette = "Set3")+
  labs(x = "Host Order 2", y = "viral similarity")+
  main_theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# temp_trefle_ntw = uni_ntw_trefle
# temp_trefle_ntw = temp_trefle_ntw[df$ID, df$ID]
# colnames(temp_trefle_ntw) = df$Order
# temp_trefle_ntw = temp_trefle_ntw%*%temp_trefle_ntw#%*%temp_trefle_ntw%*%temp_trefle_ntw%*%
#   #temp_trefle_ntw%*%temp_trefle_ntw%*%temp_trefle_ntw%*%temp_trefle_ntw
# temp_trefle_ntw = temp_trefle_ntw[1:length(unique(trefle$host)), 1:length(unique(trefle$host))]
# 
# temp = colSums(temp_trefle_ntw)
# tapply(temp, names(temp), mean)

##### direct viral sharing (A²) vs communinicability viral similarity (delta G)
ntw_trefle_sqrt = uni_ntw_trefle%*%uni_ntw_trefle
ntw_trefle_sqrt_host = ntw_trefle_sqrt[rownames(ntw_trefle_sqrt)%in% unique(trefle$host),
                                       colnames(ntw_trefle_sqrt)%in% unique(trefle$host)]
df_direct_viral_sharing = melt(ntw_trefle_sqrt_host, value.name = "dir_vrl_shr")
df_direct_viral_sharing$ID_hosts = paste(df_direct_viral_sharing$Var1, df_direct_viral_sharing$Var2)
df_direct_viral_sharing = df_direct_viral_sharing[,!(colnames(df_direct_viral_sharing) %in% c("Var1", "Var2"))]
df_clustering_host = df_clustering_host%>%
  full_join(df_direct_viral_sharing, by = c("ID_hosts"))
g3 = ggplot()+
  #geom_bin_2d(data = df_clustering_host, aes(dir_vrl_shr, clustering), bins = 50)+
  geom_smooth(data = df_clustering_host, aes(dir_vrl_shr, clustering), method = "lm")+
  geom_point(data = df_clustering_host,aes(dir_vrl_shr, clustering))+
  facet_wrap(~HostOrder2)+
  geom_hline(yintercept  = 0, col = "red", cex = 1.5, linetype =2)+
  labs(x = "Viral sharing count", y =  "Viral (spectrum) similaity")+
  main_theme
g4 = ggplot(df_clustering_host)+
  geom_bin_2d(aes(dir_vrl_shr, G), bins = 50)+
  geom_smooth(aes(dir_vrl_shr, G), method = "lm")+
  labs(x = "Viral sharing count", y =  "Communicability")+
  main_theme
ggarrange(g3, g4)

viral_sharing_mean = colMeans(ntw_trefle_sqrt_host)
viral_similarity_mean = colMeans(G_cluster_host)
df_viral_means  = data.frame(viral_sharing_mean = viral_sharing_mean,
                             viral_similarity_mean = viral_similarity_mean,
                             sp = attr(viral_sharing_mean,"names"))

df_viral_means$HostOrder = "_"
for(h in unique(df_viral_means$sp)){
  HO = clover$HostOrder[which(clover$Host == h)[1]]
  host_indexs = which(df_viral_means$sp == h)
  if (HO %in% main_h_order$HostOrder){
    df_viral_means$HostOrder[host_indexs] = as.character(HO)
  }else{
    df_viral_means$HostOrder[host_indexs] = "Other"
  }
}
str(df_viral_means)
ggplot(df_viral_means)+
  geom_point(aes(viral_sharing_mean, viral_similarity_mean))+
  geom_smooth(aes(viral_sharing_mean, viral_similarity_mean), method = "lm")+
  facet_wrap(~HostOrder)+
  main_theme
str(df_clustering_host)
df_viral_median_detail = df_clustering_host%>%
  group_by(Var1, HostOrder2)%>%
  summarise(clustering_median = median(clustering),
            G_median = median(G),
            dir_vrl_shr_median = median(dir_vrl_shr), 
            phylodist_median = median(phylodist),
            HostOrder1 = first(HostOrder1))
str(df_viral_median_detail)
ggplot(df_viral_median_detail)+
  #geom_point(aes(dir_vrl_shr_median, clustering_median , col = HostOrder2))+
  geom_smooth(aes(dir_vrl_shr_median, clustering_median, col = HostOrder2 ), se = F, method = "lm")+
  scale_colour_brewer(palette = "Set3")+
  facet_wrap(~HostOrder1)+
  labs(x = "Viral sharing count", y =  "Viral (spectrum) similaity")+
  main_theme
df_clustering_host
ggplot(df_clustering_host)+
  #geom_point(aes(dir_vrl_shr, clustering_mean , col = HostOrder2))+
  geom_smooth(aes(dir_vrl_shr, clustering, col = HostOrder2 ), se = F, method = "lm")+
  scale_colour_brewer(palette = "Set3")+
  facet_wrap(~HostOrder1)+
  main_theme

### Same for virus (host sharing)

ntw_trefle_sqrt_host = ntw_trefle_sqrt[rownames(ntw_trefle_sqrt)%in% unique(trefle$virus),
                                       colnames(ntw_trefle_sqrt)%in% unique(trefle$virus)]
df_direct_host_sharing = melt(ntw_trefle_sqrt_host, value.name = "dir_hst_shr")
df_direct_host_sharing$ID_virus = paste(df_direct_host_sharing$Var1, df_direct_host_sharing$Var2)
df_direct_host_sharing = df_direct_host_sharing[,!(colnames(df_direct_host_sharing) %in% c("Var1", "Var2"))]

df_clustering_virus = melt(G_cluster[rownames(G_cluster)%in% unique(trefle$virus),
                          colnames(G_cluster)%in% unique(trefle$virus)],value.name = "clustering")
df_clustering_virus$ID_virus = paste(df_clustering_virus$Var1, df_clustering_virus$Var2)
str(df_clustering_virus)
df_clustering_virus = df_clustering_virus%>%
  full_join(df_direct_host_sharing, by = c("ID_virus"))

main_v_order = clover%>%
  count(VirusOrder)%>%
  arrange(desc(n))%>%
  slice(1:10)
str(df_clustering_virus)
df_clustering_virus$VirusOrder1 = "_"
for(v in unique(df_clustering_virus$Var1)){
  VO = clover$VirusOrder[which(clover$Virus == v)[1]]
  virus_indexs = which(df_clustering_virus$Var1 == v)
  if (VO %in% main_v_order$VirusOrder){
    df_clustering_virus$VirusOrder1[virus_indexs] = as.character(VO)
  }else{
    df_clustering_virus$VirusOrder1[virus_indexs] = "Other"
  }
  
}
df_clustering_virus$VirusOrder2 = "_"
for(v in unique(df_clustering_virus$Var2)){
  VO = clover$VirusOrder[which(clover$Virus == v)[1]]
  virus_indexs = which(df_clustering_virus$Var2 == v)
  if (VO %in% main_v_order$VirusOrder){
    df_clustering_virus$VirusOrder2[virus_indexs] = as.character(VO)
  }else{
    df_clustering_virus$VirusOrder2[virus_indexs] = "Other"
  }
  
}
str(df_clustering_virus)
ggplot(df_clustering_virus)+
  geom_point(aes(dir_hst_shr, clustering))+
  geom_smooth(aes(dir_hst_shr, clustering), method = "lm")+
  scale_colour_brewer(palette = "Set3")+
  #facet_wrap(~HostOrder1)+
  labs(x = "Host sharing count", y =  "Host (spectrum) similaity")+
  main_theme
ggplot(df_clustering_virus)+
  geom_point(aes(dir_hst_shr, clustering, col = VirusOrder2))+
  geom_smooth(aes(dir_hst_shr, clustering, col = VirusOrder2), method = "lm")+
  scale_colour_brewer(palette = "Set3")+
  facet_wrap(~VirusOrder1)+
  labs(x = "Host sharing count", y =  "Host (spectrum) similaity")+
  main_theme
##### geo overlap data 

load("data/Finaldf.Rdata")
str(FinalHostMatrix)
FinalHostMatrix$Sp = sub("_", " ", FinalHostMatrix$Sp)
FinalHostMatrix$Sp2 = sub("_", " ", FinalHostMatrix$Sp2)
G_cluster_nrmlz = normalized_2(G_cluster)
df_G_host_temp = melt(G_cluster_nrmlz[rownames(G_cluster) %in% unique(FinalHostMatrix$Sp),
                               colnames(G_cluster) %in% unique(FinalHostMatrix$Sp2)],value.name = "clustering")
df_G_host_temp$ID_hosts = paste(df_G_host_temp$Var1, df_G_host_temp$Var2)
df_G_host_temp = df_G_host_temp[,!(colnames(df_G_host_temp) %in% c("Var1", "Var2"))]

FinalHostMatrix$ID_hosts = paste(FinalHostMatrix$Sp, FinalHostMatrix$Sp2)
FinalHostMatrix = FinalHostMatrix%>%
  full_join(df_G_host_temp, by = "ID_hosts")

ggplot(FinalHostMatrix)+
  geom_bin_2d(aes(Phylo, clustering),bins = 50)+
  geom_smooth(aes(Phylo, clustering, col = hOrder), method = "glm")+
  facet_wrap(~hOrder)+
  labs(x = "Phylogenetic similarity", y = "viral similarity")+
  main_theme
FinalHostMatrix%>%
  filter(hOrder == "PRIMATES")%>%
ggplot()+
  geom_point(aes(Phylo, clustering,col = hOrder.Sp2),bins = 50)+
  geom_smooth(aes(Phylo, clustering))+
  labs(x = "Phylogenetic similarity", y = "viral similarity")+
  main_theme
ggplot(FinalHostMatrix)+
  geom_bin_2d(aes(clustering, Space),bins = 50)+
  geom_smooth(aes(clustering, Space, col = hOrder), method = "lm")+
  facet_wrap(~hOrder)

ggplot(FinalHostMatrix)+
  geom_bin_2d(aes(clustering, Space),bins = 50)+
  geom_smooth(aes(clustering, Space, col = hOrder))+
  facet_wrap(~hOrder)
