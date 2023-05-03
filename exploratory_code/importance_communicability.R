library(tidyverse)
library(igraph)
setwd("~/Fac/Master_Rennes/stage1/edges")
source("exploratory_code/functions.R")

clover = read.csv("data/clover.csv", stringsAsFactors = T)
trefle = read.csv("data/trefle.csv", stringsAsFactors = T)

##### 1000 interaction for testing /!\to be removed/!\
set.seed(1000)
ID = sample(nrow(clover), 1000)
clover = clover[ID, ]

trefle = trefle[trefle$host %in% clover$Host & trefle$virus %in% clover$Virus ,]
#####

ID_trefle = which( paste(trefle$virus, trefle$host) %in% 
                     paste(clover$Virus, clover$Host))
shared_asso = trefle[ID_trefle, ]
str(shared_asso)
uni_ntw_clover = matrix.associations.uni(Virus = clover$Virus, Host = clover$Host)
uni_ntw_trefle = matrix.associations.uni(Virus = trefle$virus, Host = trefle$host)

par(mfrow =c(1,2))
grh_clover = graph_from_adjacency_matrix(uni_ntw_clover,mode = "undirected")
plot(grh_clover,edge.arrow.size=.2,vertex.label=NA, vertex.size=.01)
grh_trefle= graph_from_adjacency_matrix(uni_ntw_trefle,mode = "undirected")
plot(grh_trefle,edge.arrow.size=.2,vertex.label=NA, vertex.size=.01)

G_trefle = communicability(uni_ntw_trefle)
G_clover = communicability(uni_ntw_clover)
str(G_clover)
importance_trefle = vector(length = nrow(shared_asso))             
importance_clover = vector(length = nrow(shared_asso))     
for(n in 1: nrow(shared_asso)){
  ID_virus = which(row.names(uni_ntw_clover) == shared_asso$virus[n])  
  ID_host = which(colnames(uni_ntw_clover) == shared_asso$host[n]) 
  
  uni_ntw_trefle[ID_virus,ID_host] = 0 # change the interaction (perturbation == zeta)
  uni_ntw_trefle[ID_host,ID_virus] = 0 
  uni_ntw_clover[ID_virus,ID_host] = 0
  uni_ntw_clover[ID_host,ID_virus] = 0
  G_zeta_trefle = communicability(uni_ntw_trefle)
  G_zeta_clover = communicability(uni_ntw_clover)
  uni_ntw_trefle[ID_virus,ID_host] = 1  # remove the change
  uni_ntw_trefle[ID_host,ID_virus] = 1 
  uni_ntw_clover[ID_virus,ID_host] = 1
  uni_ntw_clover[ID_host,ID_virus] = 1
  
  G_delta_trefle = G_trefle-G_zeta_trefle
  G_delta_clover = G_clover-G_zeta_clover[ID_host,ID_virus]
  ###### /!\ might need a normalzation step here to compare both G_delta
  G_delta_trefle = normalized(G_delta_trefle)
  G_delta_clover = normalized(G_delta_clover)
  importance_trefle[n] = G_delta_trefle[ID_virus,ID_host]
  importance_clover[n] = G_delta_clover[ID_virus,ID_host]
  #share_eigen_T_C = c(E_temp_trefle = E_temp_trefle, E_temp_clover = E_temp_clover)
  
}
importance_df = data.frame(importance_trefle,importance_clover,shared_asso)
str(importance_df)
write.csv(importance_df, "importance_df")
ggplot(importance_df)+
  geom_point(aes(importance_clover,importance_trefle))


