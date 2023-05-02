library(tidyverse)

source("exploratory_code/functions.R")

clover = read.csv("data/clover.csv", stringsAsFactors = T)
trefle = read.csv("data/trefle.csv", stringsAsFactors = T)

##### 1000 interaction for testing /!\to be removed/!\
set.seed(1000)
ID = sample(nrow(clover), 1000)
clover_subset = clover[ID, ]
#####


uni_ntw_clover = matrix.associations.uni(Virus = clover$Virus, Host = clover$Host)
uni_ntw_trefle = matrix.associations.uni(Virus = trefle$virus, Host = trefle$host)
share_eigen_T_C = vector("list", length = nrow(shared_asso))
for(n in 1: nrow(shared_asso)){
  ID_virus = which(row.names(ntw_trefle) == shared_asso$virus[n])  
  ID_host = which(colnames(ntw_trefle) == shared_asso$host[n]) 
  
  uni_ntw_trefle[ID_virus,ID_host] = 0 # change the interaction
  uni_ntw_trefle[ID_host,ID_virus] = 0 
  uni_ntw_clover[ID_virus,ID_host] = 0
  uni_ntw_clover[ID_host,ID_virus] = 0
  E_temp_trefle = eigen(uni_ntw_trefle)
  E_temp_clover = eigen(uni_ntw_clover)
  uni_ntw_trefle[ID_virus,ID_host] = 1  # remove the change
  uni_ntw_trefle[ID_host,ID_virus] = 1 
  uni_ntw_clover[ID_virus,ID_host] = 1
  uni_ntw_clover[ID_host,ID_virus] = 1
  share_eigen_T_C[[n]] = c(E_temp_trefle = E_temp_trefle, E_temp_clover = E_temp_clover)
  
}

str(share_eigen_T_C[[n]])
str(G_trefle_eta)
str(spectra)
n =1
for(n in 1:length(share_eigen_T_C)){
  G_trefle_eta = communicability(A=uni_ntw_clover, spectra = 
                                   list(values= share_eigen_T_C[[n]]$E_temp_trefle.values ,
                                     vectors = share_eigen_T_C[[n]]$E_temp_trefle.vectors))
  G_clover = communicability(A= uni_ntw_clovers, pectra = share_eigen_T_C[[n]]$E_temp_trefle)
}




