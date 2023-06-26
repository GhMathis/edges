library(tidyverse)
library(doParallel) # for parallel processing
library(parallel) # to detect core
source("R/functions.R")

clover = read.csv("data/clover.csv", stringsAsFactors = T)
trefle = read.csv("data/trefle.csv", stringsAsFactors = T)

##### 1000 interaction for testing /!\to be removed/!\
# set.seed(1000)
# ID = sample(nrow(clover), 1000)
# clover = clover[ID, ]
# 
# trefle = trefle[trefle$host %in% clover$Host & trefle$virus %in% clover$Virus ,]
#####

ID_trefle = which( paste(trefle$virus, trefle$host) %in% 
                     paste(clover$Virus, clover$Host))
shared_asso = trefle[ID_trefle, ]

uni_ntw_clover = matrix.associations.uni(Virus = clover$Virus, Host = clover$Host)
uni_ntw_trefle = matrix.associations.uni(Virus = trefle$virus, Host = trefle$host)
all(unique(colnames(uni_ntw_trefle)) %in% unique(colnames(uni_ntw_clover)))
G_trefle = communicability(uni_ntw_trefle)
G_clover = communicability(uni_ntw_clover)
ntw_base = list(uni_ntw_clover = uni_ntw_clover, uni_ntw_trefle = uni_ntw_trefle,
                G_trefle = G_trefle, G_clover = G_clover)

communi.func<- function(virus,host,row_ID,arg){
  
  print(row_ID)
  
  with(arg,{
    
    ID_virus = which(row.names(uni_ntw_clover) == virus)  
    ID_host = which(colnames(uni_ntw_clover) == host) 
    ID_virus_t = which(row.names(uni_ntw_trefle) == virus)  
    ID_host_t = which(colnames(uni_ntw_trefle) == host)
    uni_ntw_trefle[ID_virus_t,ID_host_t] = 0 # change the interaction (perturbation == zeta)
    uni_ntw_trefle[ID_host_t,ID_virus_t] = 0 
    uni_ntw_clover[ID_virus,ID_host] = 0
    uni_ntw_clover[ID_host,ID_virus] = 0
    
    G_zeta_trefle = communicability(uni_ntw_trefle)
    
    G_zeta_clover = communicability(uni_ntw_clover)
    
    uni_ntw_trefle[ID_virus_t,ID_host_t] = 1  # remove the change
    uni_ntw_trefle[ID_host_t,ID_virus_t] = 1 
    uni_ntw_clover[ID_virus,ID_host] = 1
    uni_ntw_clover[ID_host,ID_virus] = 1
    
    ##
    G_trefle_nrmlz = normalized(G_trefle)
    G_zeta_trefle_nrmlz = normalized(G_zeta_trefle)
    G_clover_nrmlz = normalized(G_clover)
    G_zeta_clover_nrmlz = normalized(G_zeta_clover)
    
    G_delta_trefle_nrmlz = G_trefle_nrmlz-G_zeta_trefle_nrmlz
    G_delta_clover_nrmlz = G_clover_nrmlz-G_zeta_clover_nrmlz

    ##
    importance_trefle_nrmlz = G_delta_trefle_nrmlz[ID_virus_t,ID_host_t]
    importance_clover_nrmlz = G_delta_clover_nrmlz[ID_virus,ID_host]
    
    G_pq_trefle_nrmlz = G_trefle_nrmlz[ID_virus_t,ID_host_t]
    G_pq_clover_nrmlz = G_clover_nrmlz[ID_virus,ID_host]
    
    
    df =data.frame(virus = virus, host = host,
                   importance_trefle_nrmlz = importance_trefle_nrmlz,
                   importance_clover_nrmlz = importance_clover_nrmlz,
                   G_pq_trefle_nrmlz = G_pq_trefle_nrmlz,
                   G_pq_clover_nrmlz = G_pq_clover_nrmlz)

    row.names(df) = row_ID
    
    return(df)
    
  })
  
}
n_iteration = nrow(shared_asso)

ncores =detectCores()
print(ncores)
registerDoParallel(cores=(ncores-1))# Shows the number of Parallel Workers to be used

getDoParWorkers()# number of actual workers


importance_df = foreach(virus = shared_asso$virus, host = shared_asso$host, row_ID = row.names(shared_asso), arg = lapply(1:n_iteration, function (x) ntw_base),.combine = rbind, .verbose = T) %dopar%{communi.func(virus = virus, host = host, row_ID = row_ID, arg=arg)}

write.csv(importance_df, "output/importance_df.csv")


