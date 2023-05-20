library(tidyverse)
library(lattice)
library(igraph)
library(doParallel)
library(FactoMineR)
source("R/functions.R")

trefle = read.csv("data/trefle.csv", stringsAsFactors = T)
clover = read.csv("data/clover.csv", stringsAsFactors = T)
#####
set.seed(1000)
ID = sample(nrow(trefle), 1000)
trefle = trefle[ID, ]
#####
trefle$HostOrder = "_"
trefle$VirusOrder = "_"
for(n in 1:nrow(clover)){
  h = which(trefle$host == clover$Host[n])
  v = which(trefle$virus == clover$Virus[n])
  
  trefle$HostOrder[h] = as.character(clover$HostOrder[n])
  trefle$VirusOrder[v] = as.character(clover$VirusOrder[n])
}
str(trefle)
trefle_ntw = matrix.associations(Virus = trefle$virus, Host = trefle$host)
colSums(trefle_ntw)

##### unique pair of host-virus that is excusive for both in trefle
which(trefle_ntw[,colnames(trefle_ntw) =="Necromys amoenus"]==1)
which(trefle_ntw[rownames(trefle_ntw) =="Pampa virus"]==1,)
colnames(trefle_ntw[,664:665])
trefle_ntw = trefle_ntw[,-which(trefle_ntw[,colnames(trefle_ntw) =="Necromys amoenus"]==1)]

#####

host_sharing = t(trefle_ntw)%*%trefle_ntw
D = diag(diag(host_sharing)**-1/2)

diag(host_sharing) = 0
host_sharing_norm = D %*% host_sharing %*% D
host_sharing_norm[1:10,1:10]
spectra = eigen(host_sharing_norm)

levelplot(host_sharing)

levelplot(exp(host_sharing_norm))
str(spectra)

G_clusering = function(value,vector,host_names){

    temp = vector%*%t(vector)*exp(value)
    colnames(temp) = host_names
    rownames(temp) = host_names
    return(temp)
}

ncores =detectCores()
print(ncores)
registerDoParallel(cores=4)
getDoParWorkers()# Shows the number of Parallel Workers to be used

##### compute the 2 to 5 dim of communicability
clustering_list = foreach(value= spectra$values[2:5], vector= spectra$vectors[,2:5],
                          host_names= sapply(host_sharing, function(x) colnames(host_sharing)) , .verbose = T) %dopar%
  {G_clusering(value = value, vector = vector,host_names = host_names)}

str(clustering_list)

levelplot(clustering_list[[1]])
##### array to data frame

n_element = (nrow(clustering_list[[1]])*(nrow(clustering_list[[1]])-1))/2
for (ndim in 1:4){
  i= 0
  name1 = vector(mode = "numeric",n_element)
  name2 =  vector(mode = "numeric",n_element)
  value =  vector(mode = "numeric",n_element)
  
  for(r in 2:nrow(clustering_list[[ndim]])-1){
    for(c in (r+1):ncol(clustering_list[[ndim]])){
      i = i+1
      value[i] = clustering_list[[ndim]][r,c]
      name1[i] = colnames(clustering_list[[ndim]])[c]
      name2[i] = rownames(clustering_list[[ndim]])[r]
    }
  }
  if(ndim == 1){
    clusering_df = data.frame(name1,name2,value_dim2 = value)
  }else{
    clusering_df = cbind(clusering_df, value)
    names(clusering_df)[2+ndim] = paste("value_dim", ndim+1, sep="")
  }
}
clusering_df$HostOrder1= "_"
clusering_df$HostOrder2= "_"
clusering_df$name1 = as.factor(clusering_df$name1)
clusering_df$name2 = as.factor(clusering_df$name2)
for(h in as.character(levels(clusering_df$name1))){
  HO = clover$HostOrder[which(clover$Host == h)[1]]
  clusering_df$HostOrder1[clusering_df$name1 == h]= as.character(HO)
}
for(h in as.character(levels(clusering_df$name2))){
  HO = clover$HostOrder[which(clover$Host == h)[1]]
  clusering_df$HostOrder2[clusering_df$name2 == h]= as.character(HO)
}
str(clusering_df)
clusering_df$hostcompaire = clusering_df$HostOrder1 == clusering_df$HostOrder2

ggplot(clusering_df)+
  geom_point(aes(value_dim2, value_dim3, col = HostOrder1), cex =2)+
  geom_point(aes(value_dim2, value_dim3, col = HostOrder2),cex = 1.5, alpha =0.7)
ggplot(clusering_df)+
  geom_boxplot(aes(hostcompaire, value_dim2))
ggplot(clusering_df)+
  geom_point(aes(value_dim4, value_dim5,  col = HostOrder1))
str(clusering_df)

pca_dims = PCA(clusering_df[,3:6])
fviz_pca_biplotpca_dims()