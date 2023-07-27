library(tidyverse)
library(lattice)
library(igraph)
library(corrplot)
library(colorRamps)
library(RColorBrewer)
library(scico)
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

G_cluster_clover = clustering.func(uni_ntw_clover)

df = data.frame(ID = 1 : ncol(G_cluster),
                sp = colnames(G_cluster),
                partite = c(rep("host",length(unique( trefle$host))),
                            rep("virus",length(unique( trefle$virus)))
                )
)

df$Order= "_"
df$Family= "_"
# attribue host order
for(i in 1:nrow(df)){
  if(df$partite[i] == "host"){
    HO = clover$HostOrder[which(clover$Host == df$sp[i])[1]]
    df$Order[i] = as.character(HO)
    HF = clover$HostFamily[which(clover$Host == df$sp[i])[1]]
    df$Family[i] = as.character(HF)
  }else{
    VO = clover$VirusOrder[which(clover$Virus == df$sp[i])[1]]
    df$Order[i] = as.character(VO)
    VF = clover$VirusFamily[which(clover$Virus == df$sp[i])[1]]
    df$Family[i] = as.character(VF)
  }
}

df = df %>%
  arrange(partite, Order, ID)
G_cluster = normalized_2(G_cluster)
G_cluster_desc = G_cluster[df$ID, df$ID]
label_plot = df$Order
label_plot[duplicated(df$Order)]=""

make.color.scale = function(mtx, white = 0.01){
  brk <- do.breaks(c(-max(abs(mtx)), max(abs(mtx))), 8)
  pal = scico(length(brk)-1, palette = "vik")
  brk = append(brk, white, 5)
  brk = append(brk, -white, 4)
  
  pal = append(pal, "#FFFFFF", 4)
  pal = append(pal, "#FFFFFF", 4)
  
  return(list(brk = brk, pal = pal))
  
}

color_scale1 = make.color.scale(G_cluster)

levelplot(G_cluster_desc,
          scales = list(labels = label_plot, x=list(rot=90)),
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal)

#### cluster host

G_cluster_host = G_cluster[rownames(G_cluster)%in% unique(trefle$host),
                           colnames(G_cluster)%in% unique(trefle$host)]

df_host = data.frame(val = colSums(G_cluster_host), ID = 1 : ncol(G_cluster_host),
                     sp = colnames(G_cluster_host))

df_host$Order= "_"
# attribue host order
for(i in 1:nrow(df_host)){
  HO = clover$HostOrder[which(clover$Host == df_host$sp[i])[1]]
  df_host$Order[i] = as.character(HO)
}
df_host = df_host %>%
  arrange(Order)
label_plot_host = df_host$Order
label_plot_host[duplicated(df_host$Order)]=""
G_cluster_host = G_cluster_host[df_host$ID,df_host$ID]
levelplot(G_cluster_host,
          scales = list(labels =label_plot_host, x=list(rot=90)),
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal)

## recap matrix 
recap_host  = matrix(NA, ncol = length(unique(df_host$Order)), nrow =length(unique(df_host$Order)))
colnames(recap_host) = unique(df_host$Order)
rownames(recap_host) = unique(df_host$Order)
for(HO in unique(df_host$Order)){
  ID = df_host%>%
    filter(Order == HO)%>%
    pull(sp)
  for(HO2 in unique(df_host$Order)){
    ID2 = df_host%>%
      filter(Order == HO2)%>%
      pull(sp)
    recap_host[rownames(recap_host) == HO2,
               colnames(recap_host) == HO] = median(G_cluster_host[rownames(G_cluster_host) %in% ID2,
                                                                   colnames(G_cluster_host) %in% ID])
  }
}
recap_host = recap_host[order(colnames(recap_host)),order(colnames(recap_host))]

color_scale2 = make.color.scale(recap_host, white = 0.001)

levelplot(recap_host,  scales = list(x=list(rot=45)),
          at = color_scale2$brk,
          colorkey=list(at = color_scale2$brk, col = color_scale2$pal),
          col.regions = color_scale2$pal)

main_Order = df_host%>%
  count(Order)%>%
  arrange(desc(n))%>%
  slice(1:10)%>%
  pull(Order)
levelplot(recap_host[colnames(recap_host) %in% main_Order, rownames(recap_host) %in% main_Order],
          scales = list(x=list(rot=45)),
          at = color_scale2$brk,
          colorkey=list(at = color_scale2$brk, col = color_scale2$pal),
          col.regions = color_scale2$pal)
#### cluster virus


G_cluster_virus = G_cluster[rownames(G_cluster)%in% unique(trefle$virus),
                            colnames(G_cluster)%in% unique(trefle$virus)]

df_virus = data.frame(val = colSums(G_cluster_virus), ID = 1 : ncol(G_cluster_virus),
                     sp = colnames(G_cluster_virus))

df_virus$Order= "_"
df_virus$Family= "_"
# attribue virus order
str(df_virus)
for(i in 1:nrow(df_virus)){
  VO = clover$VirusOrder[which(clover$Virus == df_virus$sp[i])[1]]
  VF = clover$VirusFamily[which(clover$Virus == df_virus$sp[i])[1]]
  df_virus$Order[i] = as.character(VO)
  df_virus$Family[i] = as.character(VF)
}
str(df_virus)
df_virus = df_virus %>%
  arrange(Order)
label_plot_virus = df_virus$Order
label_plot_virus[duplicated(df_virus$Order)]=""
G_cluster_virus_plot = G_cluster_virus[df_virus$ID,df_virus$ID]



levelplot(G_cluster_virus_plot, scales = list(labels =label_plot_virus,
                                         x=list(rot=90)),
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal)

## recap matrix 
recap_virus  = matrix(NA, ncol = length(unique(df_virus$Order)), nrow =length(unique(df_virus$Order)))
colnames(recap_virus) = unique(df_virus$Order)
rownames(recap_virus) = unique(df_virus$Order)
for(VO in unique(df_virus$Order)){
  ID = df_virus%>%
    filter(Order == VO)%>%
    pull(sp)
  for(VO2 in unique(df_virus$Order)){

    ID2 = df_virus%>%
      filter(Order == VO2)%>%
      pull(sp)
    recap_virus[rownames(recap_virus) == VO2,
               colnames(recap_virus) == VO] = median(G_cluster_virus[rownames(G_cluster_virus) %in% ID2,
                                                                   colnames(G_cluster_virus) %in% ID])
  }
}
recap_virus = recap_virus[order(rownames(recap_virus)),order(colnames(recap_virus))]


levelplot(recap_virus,  scales = list(x=list(rot=45)),
          at = color_scale2$brk,
          colorkey=list(at = color_scale2$brk, col = color_scale2$pal),
          col.regions = color_scale2$pal)

main_Order_v = df_virus%>%
  count(Order)%>%
  arrange(desc(n))%>%
  slice(1:17)%>%
  pull(Order)
recap_main_virus = recap_virus[rownames(recap_virus) %in% main_Order_v,
            colnames(recap_virus) %in% main_Order_v]
levelplot(recap_main_virus,
          scales = list(x=list(rot=45)),
          at = color_scale2$brk,
          colorkey=list(at = color_scale2$brk, col = color_scale2$pal),
          col.regions = color_scale2$pal)

#### cluster host virus


G_cluster_hv = G_cluster[rownames(G_cluster)%in% unique(trefle$virus),
                            colnames(G_cluster)%in% unique(trefle$host)]

G_cluster_hv = G_cluster_hv[df_virus$ID,df_host$ID]

G_cluster_hv_plot = G_cluster_hv
colnames(G_cluster_hv_plot) =label_plot_host
rownames(G_cluster_hv_plot) =label_plot_virus

levelplot(G_cluster_hv_plot,
          scales = list(x=list(rot=90)),
          xlab = "virus order",
          ylab = "host_order",
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal)

## recap matrix 
recap_hv  = matrix(NA, nrow = length(unique(df_virus$Order)), ncol =length(unique(df_host$Order)))
rownames(recap_hv) = unique(df_virus$Order)
colnames(recap_hv) = unique(df_host$Order)
for(HO in unique(df_host$Order)){
  ID = df_host%>%
    filter(Order == HO)%>%
    pull(sp)
  for(VO in unique(df_virus$Order)){
    ID2 = df_virus%>%
      filter(Order == VO)%>%
      pull(sp)

    recap_hv[rownames(recap_hv) == VO,
                colnames(recap_hv) == HO] = median(G_cluster_hv[rownames(G_cluster_hv) %in% ID2,
                                                                      colnames(G_cluster_hv) %in% ID])
  }
}
recap_hv = recap_hv[order(rownames(recap_hv)), order(colnames(recap_hv))]


color_scale3 = make.color.scale(recap_hv, white = 0.001)
levelplot(recap_hv,  scales = list(x=list(rot=45)),
          xlab = "virus order",
          ylab = "host_order",
          at = color_scale3$brk,
          colorkey=list(at = color_scale3$brk, col = color_scale3$pal),
          col.regions = color_scale3$pal)


levelplot(recap_hv[rownames(recap_hv) %in% main_Order_v,colnames(recap_hv) %in% main_Order],
          xlab = "virus order",
          ylab = "host_order",
          at = color_scale3$brk,
          colorkey=list(at = color_scale3$brk, col = color_scale3$pal),
          col.regions = color_scale3$pal)

### compaire to random
set.seed(1000)
rand_edges_ntw = matrix.associations.uni(sample(trefle$virus),trefle$host)
G_cluster_rand = clustering.func(rand_edges_ntw)
str(G_cluster_rand)
G_cluster_rand = normalized_2(G_cluster_rand)

color_scale4 = make.color.scale(G_cluster_rand)

levelplot(G_cluster_rand[df$ID, df$ID],  scales = list(labels = label_plot, x=list(rot=45)),
          xlab = "virus order",
          ylab = "host_order",
          at = color_scale4$brk,
          colorkey=list(at = color_scale4$brk, col = color_scale4$pal),
          col.regions = color_scale4$pal)


##### Comparison with Albery

PredictedNetwork_host = readRDS("data/PredictedNetwork.rds")
str(PredictedNetwork_host)
colnames(PredictedNetwork_host) = sub("_", " ", colnames(PredictedNetwork_host))
rownames(PredictedNetwork_host) = sub("_", " ", rownames(PredictedNetwork_host))

host_order_df = read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt")
host_order_df = host_order_df[,c(1,3,4)]
host_order_df$names_host = paste(host_order_df$MSW05_Genus, host_order_df$MSW05_Species,sep=" ")
host_order_df = host_order_df %>%
  arrange(MSW05_Order)
PredictedNetwork_host  =as.matrix(PredictedNetwork_host)
host_order_df$ID = NA 
for(h in colnames(PredictedNetwork_host)){
  host_order_df$ID[host_order_df$names_host == h] = which(colnames(PredictedNetwork_host)==h)
}
str(host_order_df)
host_order_df = na.omit(host_order_df)

label_plot_host2 = host_order_df$MSW05_Order
label_plot_host2[duplicated(host_order_df$MSW05_Order)]=""

levelplot(PredictedNetwork_host[host_order_df$ID,host_order_df$ID],
          scales = list(labels = label_plot_host2, x=list(rot=90)))


recap_sharing  =matrix(NA, ncol = length(unique(host_order_df$MSW05_Order)), nrow =length(unique(host_order_df$MSW05_Order)))
colnames(recap_sharing) = unique(host_order_df$MSW05_Order)
rownames(recap_sharing) = unique(host_order_df$MSW05_Order)
for(HO in unique(host_order_df$MSW05_Order)){
  for(HO2 in unique(host_order_df$MSW05_Order)){
    ID = host_order_df%>%
      filter(MSW05_Order == HO)%>%
      pull(ID)
    ID2 = host_order_df%>%
      filter(MSW05_Order == HO2)%>%
      pull(ID)
    recap_sharing[colnames(recap_sharing) == HO,
           rownames(recap_sharing) == HO2] = median(PredictedNetwork_host[ID,ID2])
  }
}

levelplot(recap_sharing,scales = list(x=list(rot=45)))

df_sharing_prob = melt(PredictedNetwork_host,value.name = "prob_sharing") 
df_clustering = melt(G_cluster_host,value.name = "clustering")
df_clustering_rand = melt(G_cluster_rand[rownames(G_cluster)%in% unique(trefle$host),
                                         colnames(G_cluster)%in% unique(trefle$host)],value.name = "clustering_rand")
str(df_clustering_rand)
df_compaire = df_clustering
df_compaire = cbind(df_compaire,clustering_rand = df_clustering_rand$clustering_rand )
df_sharing_prob$sharing_ID = paste(df_sharing_prob$Var1, df_sharing_prob$Var2)
df_sharing_prob = df_sharing_prob%>%
  select(prob_sharing,sharing_ID)
df_compaire = df_compaire%>%
  mutate(sharing_ID = paste(Var1,Var2))%>%
  full_join(df_sharing_prob, by = "sharing_ID", keep = F)

str(df_compaire)

ggplot(df_compaire)+
  geom_bin_2d(aes(clustering, prob_sharing),bins = 50)+
  geom_smooth(aes(clustering, prob_sharing))

ggplot(df_compaire)+
  geom_bin_2d(aes(clustering, clustering_rand),bins = 50)+
  geom_smooth(aes(clustering, clustering_rand))

ggplot(df_compaire)+
  geom_bin_2d(aes(clustering_rand, prob_sharing),bins = 50)+
  geom_smooth(aes(clustering_rand, prob_sharing))

###### clustering Albery Host Host sharing ntw
PredictedNetwork_cluster = clustering.func(PredictedNetwork_host)
str(PredictedNetwork_cluster)
PredictedNetwork_cluster = normalized_2(PredictedNetwork_cluster)
color_scale_cluster_Alb = make.color.scale(PredictedNetwork_cluster)
levelplot(PredictedNetwork_cluster[host_order_df$ID,host_order_df$ID],
          scales = list(labels = label_plot_host2, x=list(rot=90)),
          at = color_scale_cluster_Alb$brk,
          colorkey=list(at = color_scale_cluster_Alb$brk, col = color_scale_cluster_Alb$pal),
          col.regions = color_scale_cluster_Alb$pal)

##### Distrib

df_cluster_hv = melt(G_cluster_hv,value.name = "clustering")
trefle_ntw = matrix.associations(trefle$virus, trefle$host)
trefle_ntw = trefle_ntw[df_virus$ID, df_host$ID]
df_ntw = melt(trefle_ntw,value.name = "association")
df_cluster_hv = cbind(df_cluster_hv, association = df_ntw$association)
names(df_cluster_hv)[1:2] = c("virus", "host")
str(df_cluster_hv)
df_cluster_hv$association = as.factor(df_cluster_hv$association)
ggplot(df_cluster_hv)+
  geom_density(aes(x = clustering, y = after_stat(scaled), fill = association), alpha = 0.4)+
  labs(y="scaled density")+
  main_theme

###### Focus on primate #####
clover_primate = clover%>%
  filter(HostOrder == "Primates")
str(clover_primate)
G_cluster_primate = G_cluster[rownames(G_cluster)%in% unique(clover_primate$Host),
                           colnames(G_cluster)%in% unique(clover_primate$Host)]

str(G_cluster_primate)
levelplot(G_cluster_primate,
          scales = list(x=list(rot=90)),
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal)

df_primate = data.frame(val = colSums(G_cluster_primate), ID = 1 : ncol(G_cluster_primate),
                     sp = colnames(G_cluster_primate))

df_primate$Family= "_"
str(clover)
# attribue host order
for(i in 1:nrow(df_primate)){
  HF = clover$HostFamily[which(clover$Host == df_primate$sp[i])[1]]
  df_primate$Family[i] = as.character(HF)
}
df_primate = df_primate %>%
  arrange(Family)
label_plot_primate = df_primate$Family
label_plot_primate[duplicated(df_primate$Family)]=""
G_cluster_primate = G_cluster_primate[df_primate$ID,df_primate$ID]

levelplot(G_cluster_primate,
          scales = list(label =label_plot_primate ,x=list(rot=90)),
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal)



recap_primate  = matrix(NA, nrow = length(unique(df_primate$Family)),
                        ncol =length(unique(df_primate$Family)))
rownames(recap_primate) = unique(df_primate$Family)
colnames(recap_primate) = unique(df_primate$Family)
for(HF in unique(df_primate$Family)){
  ID = df_primate%>%
    filter(Family == HF)%>%
    pull(sp)
  for(HF2 in unique(df_primate$Family)){
    ID2 = df_primate%>%
      filter(Family == HF2)%>%
      pull(sp)
    
    recap_primate[rownames(recap_primate) == HF2,
             colnames(recap_primate) == HF] = median(G_cluster_primate[rownames(G_cluster_primate) %in% ID2,
                                                             colnames(G_cluster_primate) %in% ID])
  }
}
color_scale5 = make.color.scale(recap_primate)
levelplot(recap_primate,
          scales = list(x=list(rot=90)),
          at = color_scale5$brk,
          colorkey=list(at = color_scale5$brk, col = color_scale5$pal),
          col.regions = color_scale5$pal)

##### top virus global health security
sort(unique(df_virus$Family))
unique(df_virus$Order)

df_virus_top_prio = df_virus%>%
  filter(Family %in%c("Poxviridae","Picobirnaviridae","Flaviviridae",
                      "Pneumoviridae","Filoviridae","Coronaviridae",
                      "Paramyxoviridae", "Rhabdoviridae", "Togaviridae", "Reoviridae") | Order %in% c("Bunyavirales"))
str(df_virus_top_prio)
df_virus_top_prio = df_virus_top_prio%>%
  arrange(Family)

virus_top_prio_fam = df_virus_top_prio$Family
virus_top_prio_fam[duplicated(df_virus_top_prio$Family)]=""

G_virus_top_prio = G_cluster_virus[df_virus_top_prio$ID, df_virus_top_prio$ID]

colnames(G_virus_top_prio)
#label = virus_top_prio_fam,
str(df_virus_top_prio$sp)
str(G_virus_top_prio)
str(virus_top_prio_fam)
levelplot(G_virus_top_prio, scales = list(label = virus_top_prio_fam,
                                         x=list(rot=90)),
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal)
##### 

known_asso_sapiens = clover%>%
  filter(Host == "Homo sapiens")%>%
  select(c(Host, Virus, VirusOrder, VirusFamily))
str(known_asso_sapiens)
G_cluster_sapien =data.frame( Val = G_cluster[,colnames(G_cluster) == "Homo sapiens"],
                               sp = rownames(G_cluster), ID = 1:ncol(G_cluster))

G_cluster_sapien$partite = c(rep("host", length(unique(trefle$host))),
                             rep("virus", length(unique(trefle$virus))))

G_cluster_sapien = G_cluster_sapien%>%
  arrange(desc(Val), group_by =partite )

G_cluster_sapien = G_cluster_sapien%>%
  mutate(known_asso = sp %in% known_asso_sapiens$Virus)



G_cluster_sapien$Order= "_"
G_cluster_sapien$Family= "_"
# attribue host order
for(i in 1:nrow(G_cluster_sapien)){
  if(G_cluster_sapien$partite[i] == "host"){
    HO = clover$HostOrder[which(clover$Host == G_cluster_sapien$sp[i])[1]]
    G_cluster_sapien$Order[i] = as.character(HO)
    HF = clover$HostFamily[which(clover$Host == G_cluster_sapien$sp[i])[1]]
    G_cluster_sapien$Family[i] = as.character(HF)
  }else{
    VO = clover$VirusOrder[which(clover$Virus == G_cluster_sapien$sp[i])[1]]
    G_cluster_sapien$Order[i] = as.character(VO)
    VF = clover$VirusFamily[which(clover$Virus == G_cluster_sapien$sp[i])[1]]
    G_cluster_sapien$Family[i] = as.character(VF)
  }
}
str(G_cluster_sapien)

G_cluster_sapien_max_v = G_cluster_sapien%>%
  filter(!known_asso, partite == "virus")%>%
  arrange(desc(Val))%>%
  slice(1:100)%>%
  filter(sp %in% df_virus_top_prio$sp)%>%
  arrange(Order,Family)

G_cluster_sapien_max_h = G_cluster_sapien%>%
  filter(partite == "host")%>%
  arrange(desc(Val))%>%
  slice(1:25)
G_cluster_sapien_max_h = rbind(G_cluster_sapien_max_h,
                               G_cluster_sapien%>%
                                 filter(sp == "Homo sapiens"))

G_cluster_sapien_max_h = G_cluster_sapien_max_h%>%
  arrange(Order,Family)

G_cluster_sapien_max = rbind(G_cluster_sapien_max_v, G_cluster_sapien_max_h)
index_sp_names = G_cluster_sapien_max$sp

G_cluster_matrix_sapien = G_cluster[G_cluster_sapien_max$ID,
                                    G_cluster_sapien_max$ID]


labels_sapien_max = G_cluster_sapien_max$Family
labels_sapien_max[duplicated(labels_sapien_max)]=""
labels_sapien_max[G_cluster_sapien_max$sp == "Homo sapiens"] = "Homo sapiens"
levelplot(G_cluster_matrix_sapien, scales = list(label = labels_sapien_max,
          x=list(rot=90)),
          at = color_scale1$brk,
          colorkey=list(at = color_scale1$brk, col = color_scale1$pal),
          col.regions = color_scale1$pal)
