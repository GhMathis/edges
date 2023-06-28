library(tidyverse)
library(doParallel) # for parallel processing
library(parallel) # to detect core
source("R/functions.R")

#clover = read.csv("data/clover.csv", stringsAsFactors = T)
trefle = read.csv("data/trefle.csv", stringsAsFactors = T)
imputed = read.csv("data/imputed_associations.csv", stringsAsFactors = T)

##### 1000 interaction for testing /!\to be removed/!\
# set.seed(1000)
# ID = sample(nrow(imputed), 1000)
# imputed = imputed[ID, ]
#  
# trefle = trefle[trefle$host %in% imputed$host & trefle$virus %in% imputed$virus ,]
#####

trefle$virushost = paste(trefle$virus,trefle$host)
imputed$virushost = paste(imputed$virus,imputed$host)
imputed = imputed[names(imputed)%in% c("P","virushost" )]
trefle = trefle%>%
  full_join(imputed, by = c("virushost"), keep=F)
trefle$P[is.na(trefle$P)] = 1

uni_ntw_trefle_prob = matrix.associations.uni(Virus = trefle$virus,
                                              Host = trefle$host,
                                              prob = trefle$P)
uni_ntw_trefle = matrix.associations.uni(Virus = trefle$virus, Host = trefle$host)

all(unique(colnames(uni_ntw_trefle)) %in% unique(colnames(uni_ntw_trefle_prob)))

G_trefle = communicability(uni_ntw_trefle)
G_trefle_prob = communicability(uni_ntw_trefle_prob)

ntw_base = list(uni_ntw_trefle_prob = uni_ntw_trefle_prob, uni_ntw_trefle = uni_ntw_trefle,
                G_trefle = G_trefle, G_trefle_prob = G_trefle_prob)

communi.func<- function(virus,host,row_ID,arg){
  
  print(row_ID)
  
  with(arg,{
    
    ID_virus_p = which(row.names(uni_ntw_trefle_prob) == virus)  
    ID_host_p = which(colnames(uni_ntw_trefle_prob) == host) 
    ID_virus_t = which(row.names(uni_ntw_trefle) == virus)  
    ID_host_t = which(colnames(uni_ntw_trefle) == host)
    
    uni_ntw_trefle[ID_virus_t,ID_host_t] = 0 # change the interaction (perturbation == zeta)
    uni_ntw_trefle[ID_host_t,ID_virus_t] = 0 
    uni_ntw_trefle_prob[ID_virus_p,ID_host_p] = 0
    uni_ntw_trefle_prob[ID_host_p,ID_virus_p] = 0
    
    G_zeta_trefle = uni_ntw_trefle
    
    G_zeta_trefle_prob = communicability(uni_ntw_trefle_prob)
    
    uni_ntw_trefle[ID_virus_t,ID_host_t] = 1  # remove the change
    uni_ntw_trefle[ID_host_t,ID_virus_t] = 1 
    uni_ntw_trefle_prob[ID_virus_p,ID_host_p] = 1
    uni_ntw_trefle_prob[ID_host_p,ID_virus_p] = 1
    
    ##
    G_trefle_nrmlz = normalized(G_trefle)
    G_zeta_trefle_nrmlz = normalized(G_zeta_trefle)
    G_trefle_prob_nrmlz = normalized(G_trefle_prob)
    G_zeta_trefle_prob_nrmlz = normalized(G_zeta_trefle_prob)
    
    G_delta_trefle_nrmlz = G_trefle_nrmlz-G_zeta_trefle_nrmlz
    G_delta_trefle_prob_nrmlz = G_trefle_prob_nrmlz-G_zeta_trefle_prob_nrmlz

    ##
    importance_trefle_nrmlz = G_delta_trefle_nrmlz[ID_virus_t,ID_host_t]
    importance_trefle_prob_nrmlz = G_delta_trefle_prob_nrmlz[ID_virus_p,ID_host_p]
    
    G_pq_trefle_nrmlz = G_trefle_nrmlz[ID_virus_t,ID_host_t]
    G_pq_trefle_prob_nrmlz = G_trefle_prob_nrmlz[ID_virus_p,ID_host_p]
    
    
    df =data.frame(virus = virus, host = host,
                   #importance_trefle_nrmlz = importance_trefle_nrmlz,
                   importance_trefle_prob_nrmlz = importance_trefle_prob_nrmlz,
                   #G_pq_trefle_nrmlz = G_pq_trefle_nrmlz,
                   G_pq_trefle_prob_nrmlz = G_pq_trefle_prob_nrmlz)

    row.names(df) = row_ID
    
    return(df)
    
  })
  
}

print(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))

array_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n_array = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))

### find the good number of pertubation to compute in each array. 
# Use the number of perturbation to compute and the total number of slurm array
n_iteration = nrow(trefle)
n = n_iteration %/% n_array
r = n_iteration %% n_array

if(array_id == 1){
  iter = c(1:(n*(array_id)-1))
  print(head(iter,n = 1L))
  print(tail(iter,n = 1L))
}else{
  iter = c((n*(array_id-1)):(n*(array_id)-1))
  print(head(iter,n = 1L))
  print(tail(iter,n = 1L))
}
if(array_id == n_array){
  iter = c(iter,c((n*array_id):n_iteration))
  print(head(iter,n = 1L))
  print(tail(iter,n = 1L))
}


ncores =detectCores()
print(ncores)
registerDoParallel(cores=ncores-1)# Shows the number of Parallel Workers to be used
getDoParWorkers()# number of actual workers

importance_binary_proba_df = foreach(virus = trefle$virus[iter], host = trefle$host[iter], row_ID = iter, arg = lapply(1:length(iter), function (x) ntw_base),.combine = rbind, .verbose = T) %dopar%{communi.func(virus = virus, host = host, row_ID = row_ID, arg=arg)}


write.csv(importance_binary_proba_df, paste("output/importance_binary_proba_df", array_id, ".csv", sep=""))


