library(tidyverse)
library(doParallel) # for parallel processing
library(parallel) # to detect core
getwd()
#setwd("~/lustre06/project/6002172/edges")
#setwd("~/Fac/Master_Rennes/stage1/edges")
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
G_trefle = communicability(uni_ntw_trefle)
G_clover = communicability(uni_ntw_clover)
ntw_base = list(uni_ntw_clover = uni_ntw_clover, uni_ntw_trefle = uni_ntw_trefle,
                G_trefle = G_trefle, G_clover = G_clover)

communi.func<- function(virus,host,row_ID,arg){
  
  print(row_ID)
  
  with(arg,{
    
    ID_virus = which(row.names(uni_ntw_clover) == virus)  
    ID_host = which(colnames(uni_ntw_clover) == host) 
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
    G_delta_clover = G_clover-G_zeta_clover
    
    ##
    G_trefle_nrmlz = normalized(G_trefle)
    G_zeta_trefle_nrmlz = normalized(G_zeta_trefle)
    G_clover_nrmlz = normalized(G_clover)
    G_zeta_clover_nrmlz = normalized(G_zeta_clover)
    
    G_delta_trefle_nrmlz = G_trefle_nrmlz-G_zeta_trefle_nrmlz
    G_delta_clover_nrmlz = G_clover_nrmlz-G_zeta_clover_nrmlz
    
    ##
    importance_trefle = G_delta_trefle[ID_virus,ID_host]
    importance_clover = G_delta_clover[ID_virus,ID_host]
    
    f_norm_G_trefle = norm(G_zeta_trefle,"F") - norm(G_trefle,"F")
    f_norm_G_clover = norm(G_zeta_clover,"F") - norm(G_clover,"F")
    
    # G_pq_zeta_trefle = G_zeta_trefle[ID_virus,ID_host]
    # G_pq_zeta_clover = G_zeta_clover[ID_virus,ID_host]
    
    G_pq_trefle = G_trefle[ID_virus,ID_host]
    G_pq_clover = G_clover[ID_virus,ID_host]
    
    ##
    importance_trefle_nrmlz = G_delta_trefle_nrmlz[ID_virus,ID_host]
    importance_clover_nrmlz = G_delta_clover_nrmlz[ID_virus,ID_host]
    
    f_norm_G_trefle_nrmlz = norm(G_zeta_trefle_nrmlz,"F") - norm(G_trefle_nrmlz,"F")
    f_norm_G_clover_nrmlz = norm(G_zeta_clover_nrmlz,"F") - norm(G_clover_nrmlz,"F")
    
    # G_pq_zeta_trefle_nrmlz = G_zeta_trefle_nrmlz[ID_virus,ID_host]
    # G_pq_zeta_clover_nrmlz = G_zeta_clover_nrmlz[ID_virus,ID_host]
    
    G_pq_trefle_nrmlz = G_trefle_nrmlz[ID_virus,ID_host]
    G_pq_clover_nrmlz = G_clover_nrmlz[ID_virus,ID_host]
    
    
    df =data.frame(virus = virus, host = host,
                   importance_trefle = importance_trefle,importance_clover = importance_clover,
                   f_norm_G_trefle = f_norm_G_trefle, f_norm_G_clover = f_norm_G_clover,
                   G_pq_trefle = G_pq_trefle, G_pq_clover = G_pq_clover,
                   importance_trefle_nrmlz = importance_trefle_nrmlz,
                   importance_clover_nrmlz = importance_clover_nrmlz,
                   f_norm_G_trefle_nrmlz = f_norm_G_trefle_nrmlz,
                   f_norm_G_clover_nrmlz = f_norm_G_clover_nrmlz,
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


communicability_df = foreach(virus = shared_asso$virus, host = shared_asso$host, row_ID = row.names(shared_asso), arg = lapply(1:n_iteration, function (x) ntw_base),.combine = rbind, .verbose = T) %dopar%{communi.func(virus = virus, host = host, row_ID = row_ID, arg=arg)}

write.csv(communicability_df, "output/communicability_df.csv")


