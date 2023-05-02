library(tidyverse)
library(ggpubr)
library(svd)
library(rgl)
library(igraph)

source("exploratory_code/functions.R")

clover = read.csv("data/clover.csv", stringsAsFactors = T)
trefle = read.csv("data/trefle.csv", stringsAsFactors = T)

ntw_clover = matrix.associations(Virus = clover$Virus, Host = clover$Host)
ntw_trefle = matrix.associations(Virus = trefle$virus, Host = trefle$host)
##### 1000 interaction for testing /!\to be removed/!\
set.seed(1000)
ID = sample(nrow(clover), 1000)
clover = clover[ID, ]
#####

ID_trefle = which( paste(trefle$virus, trefle$host) %in% paste(clover$Virus, clover$Host))
length(ID_trefle) # lot od dipicate in clover (mutiple publications for same assocations)

shared_asso = trefle[ID_trefle, ]

I0_clover = svd.rpd(clover_ntw)
I0_trefle = svd.rpd(trefle_ntw)

O_diff_clover_L = vector(length = nrow(shared_asso))
O_diff_clover_R = vector(length = nrow(shared_asso))
O_diff_trefle_L = vector(length = nrow(shared_asso))
O_diff_trefle_R = vector(length = nrow(shared_asso))
share_svd_T_C = vector("list", length = nrow(shared_asso))

for(n in 1: nrow(shared_asso)){
  ID_virus = which(row.names(ntw_trefle) == shared_asso$virus[n])  
  ID_host = which(colnames(ntw_trefle) == shared_asso$host[n]) 
  
  ntw_trefle[ID_virus,ID_host] = 0 # change the interaction
  ntw_clover[ID_virus,ID_host] = 0
  I_temp_trefle = svd.rpd(ntw_trefle)
  I_temp_clover = svd.rpd(ntw_clover)
  ntw_trefle[ID_virus,ID_host] = 1 # remove the change
  ntw_clover[ID_virus,ID_host] = 1
  
  O_diff_clover_L[n] = norm(I_temp_clover$L, "O")- norm(I0_clover$L, "O")
  O_diff_clover_R[n] = norm(I_temp_clover$t_R, "O") - norm(I0_clover$t_R, "O")
  
  O_diff_trefle_L[n] = norm(I_temp_trefle$L, "O") - norm(I0_trefle$L, "O")
  O_diff_trefle_R[n] = norm(I_temp_trefle$t_R, "O") - norm(I0_trefle$t_R, "O")
  
  share_svd_T_C[[n]] = c(svd_trefle = I_temp_trefle, svd_clover = I_temp_clover)
}
for(n in 1: nrow(shared_asso)){
  O_diff_clover_L[n] = norm(share_svd_T_C[[n]]$svd_clover.L, "O")- norm(I0_clover$L, "O")
  O_diff_clover_R[n] = norm(share_svd_T_C[[n]]$svd_clover.t_R) - norm(I0_clover$t_R, "O")
  
  O_diff_trefle_L[n] = norm(share_svd_T_C[[n]]$svd_trefle.L, "O") - norm(I0_trefle$L, "O")
  O_diff_trefle_R[n] = norm(share_svd_T_C[[n]]$svd_trefle.t_R, "O") - norm(I0_trefle$t_R, "O")
}
all_norm = data.frame(O_diff_clover_L = O_diff_clover_L,
                      O_diff_clover_R = O_diff_clover_R,
                      O_diff_trefle_L = O_diff_trefle_L,
                      O_diff_trefle_R = O_diff_trefle_R,
                      shared_asso)
normalized = function(x) (x-min(x))/(max(x)-min(x))

all_norm = data.frame(sapply(all_norm[,1:4], normalized),shared_asso )
str(all_norm)
ggplot(all_norm)+
  geom_point(aes(O_diff_clover_L+O_diff_clover_R,O_diff_trefle_L+O_diff_trefle_R))+
  geom_abline(intercept = 0,slope = 1, col="red")
ggplot(all_norm)+
  geom_point(aes(O_diff_clover_R,O_diff_trefle_R))+
  geom_abline(intercept = 0,slope = 1, col="red")

