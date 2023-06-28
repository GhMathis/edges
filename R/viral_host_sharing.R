library(tidyverse)
library(doParallel) # for parallel processing
library(parallel) # to detect core

source("R/functions.R")

clover = read.csv("data/clover.csv", stringsAsFactors = T)
trefle = read.csv("data/trefle.csv", stringsAsFactors = T)
##### 1000 interaction for testing 
set.seed(1000)
ID = sample(nrow(clover), 1000)
clover = clover[ID, ]

trefle = trefle[trefle$host %in% clover$Host & trefle$virus %in% clover$Virus ,]

##### 
trefle$HostOrder= "_"
trefle$VirusOrder= "_"
trefle$host = as.factor(trefle$host)
trefle$virus = as.factor(trefle$virus)

# attribue host order
for(h in as.character(levels(trefle$host))){
  HO = clover$HostOrder[which(clover$Host == h)[1]]
  trefle$HostOrder[trefle$host == h]= as.character(HO)
}
# attribue virus order
for(v in as.character(levels(trefle$virus))){
  VO = clover$VirusOrder[which(clover$Virus == v)[1]]
  trefle$VirusOrder[trefle$virus == v]= as.character(VO)
}

##### compute adjacency matrix 
uni_ntw_trefle = matrix.associations.uni(Virus = trefle$virus, Host = trefle$host)
##### compute communicability matrix
G_trefle = communicability(uni_ntw_trefle)


trefle$HostOrder = as.character(trefle$HostOrder)
trefle$VirusOrder = as.character(trefle$VirusOrder)

##### list given to function to identify host and virus order + unpreturbed
# adjacency and communicability matrix

ntw_base = list(uni_ntw_trefle = uni_ntw_trefle,
                G_trefle = G_trefle,
                trefle=trefle)
virus= trefle$virus[1]
host = trefle$host[1]
communi.func<- function(virus,host,HostOrder_zeta, VirusOrder_zeta, row_ID,arg){
  
  print(row_ID)
  
  with(arg,{
    ID_virus = which(row.names(uni_ntw_trefle) == virus)  
    ID_host = which(colnames(uni_ntw_trefle) == host) 
    
    ###preturbed adjacency matrix
    uni_ntw_trefle[ID_virus,ID_host] = 0 # change the interaction (perturbation == zeta)
    uni_ntw_trefle[ID_host,ID_virus] = 0 
    
    # preturbed communicability matrix
    G_zeta_trefle = communicability(uni_ntw_trefle)
    
    # adjacency matrix
    uni_ntw_trefle[ID_virus,ID_host] = 1  # remove the change
    uni_ntw_trefle[ID_host,ID_virus] = 1
    
    
    ##### Compute all z score for all groups of Host and Virus order.
    ### HOST


    # Host in the order with the pertubation
    zeta_host_order = trefle$host[trefle$HostOrder == HostOrder_zeta]
    ID_intraOrder_h =  which(colnames(uni_ntw_trefle) %in% unique(zeta_host_order))
    
    # Host out of order with the pertubation
    zeta_host_extraorder = trefle$host[trefle$HostOrder != HostOrder_zeta]
    ID_extraOrder_h =  which(colnames(uni_ntw_trefle) %in% unique(zeta_host_extraorder))
    
    host_host_ID = c(ID_intraOrder_h,ID_extraOrder_h)
    
    # global effect of the pertubation 
    mean_dG_h = mean(G_trefle[host_host_ID, host_host_ID]-
                       G_zeta_trefle[host_host_ID, host_host_ID])
    
    sd_dG_h = sd(G_trefle[host_host_ID, host_host_ID]
                 -G_zeta_trefle[host_host_ID, host_host_ID])
    
    #compute G and G zeta means for intra and extra order
    g_zeta_intra_order_h = mean(G_zeta_trefle[ID_intraOrder_h, ID_intraOrder_h])
    g_intra_order_h = mean(G_trefle[ID_intraOrder_h, ID_intraOrder_h])
    
    
    g_zeta_extra_order_h = mean(G_zeta_trefle[ID_extraOrder_h, ID_extraOrder_h])
    g_extra_order_h = mean(G_trefle[ID_extraOrder_h, ID_extraOrder_h])
    
    # mean(G - G zeta) 
    g_d_intra_order_h = mean(G_trefle[ID_intraOrder_h,ID_intraOrder_h]-
      G_zeta_trefle[ID_intraOrder_h, ID_intraOrder_h])
    g_d_extra_order_h = mean(G_trefle[ID_extraOrder_h,ID_extraOrder_h]-
                               G_zeta_trefle[ID_extraOrder_h, ID_extraOrder_h])
    
    score_intra_h = ((g_d_intra_order_h)- 
                                 (mean_dG_h))/(sd_dG_h)
    
    score_extra_h = ((g_d_extra_order_h )-
                                 (mean_dG_h))/ (sd_dG_h)
    
    
    ### VIRUS
    # Virus in the order with the pertubation
    zeta_virus_order = trefle$virus[trefle$VirusOrder == VirusOrder_zeta]
    ID_intraOrder_v =  which(colnames(uni_ntw_trefle) %in% unique(zeta_virus_order))
    
    # Virus out of order with the pertubation
    zeta_virus_extraorder= trefle$virus[trefle$VirusOrder != VirusOrder_zeta]
    ID_extraOrder_v =  which(colnames(uni_ntw_trefle) %in% unique(zeta_virus_extraorder))
    
    virus_virus_ID = c(ID_intraOrder_v,ID_extraOrder_v)
    
    # global effect of the pertubation 
    mean_dG_v = mean(G_trefle[virus_virus_ID,virus_virus_ID]-
                          G_zeta_trefle[virus_virus_ID,virus_virus_ID])
    sd_dG_v = sd(G_zeta_trefle[virus_virus_ID,virus_virus_ID]-
                      G_trefle[virus_virus_ID,virus_virus_ID])

    
    
    #compute G and G zeta means for intra and extra order
    g_zeta_intra_order_v = mean(G_zeta_trefle[ID_intraOrder_v, ID_intraOrder_v])
    g_intra_order_v = mean(G_trefle[ID_intraOrder_v, ID_intraOrder_v])
    
    g_extra_order_v = mean(G_trefle[ID_extraOrder_v, ID_extraOrder_v])
    g_zeta_extra_order_v = mean(G_zeta_trefle[ID_extraOrder_v, ID_extraOrder_v])
    
    
    g_d_intra_order_v = mean(G_zeta_trefle[ID_intraOrder_v, ID_intraOrder_v]-
                               G_trefle[ID_intraOrder_v,ID_intraOrder_v])
    g_d_extra_order_v = mean(G_zeta_trefle[ID_extraOrder_v, ID_extraOrder_v] -
                               G_trefle[ID_extraOrder_v,ID_extraOrder_v])
    score_intra_v = ((g_d_intra_order_v)- 
                       (mean_dG_v))/(sd_dG_v)
    
    score_extra_v = ((g_d_extra_order_v )-
                       (mean_dG_v))/ (sd_dG_v)
    
    
    df = data.frame(HostOrder_zeta = HostOrder_zeta,
                    host = host,
                    score_intra_h,score_extra_h,
                    VirusOrder_zeta = VirusOrder_zeta,
                    virus = virus,
                    score_intra_v,score_extra_v
                    )
    ##
    
    
    
    
    df$X = row_ID
    
    return(df)
    
  })
  
}

##### slurm array for big computation. 
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

### test
#iter = 1:getDoParWorkers()
###

viral_host_sharing_G_df = foreach(virus = trefle$virus[iter], host = trefle$host[iter],
                                                                HostOrder_zeta = trefle$HostOrder[iter],
                                                                VirusOrder_zeta = trefle$VirusOrder[iter],
                                                                row_ID =iter,
                                                                arg = lapply(1:length(iter), function (x) ntw_base),
                                                                .combine = rbind, .verbose = T) %dopar%
                                 {communi.func(virus = virus, host = host, HostOrder_zeta = HostOrder_zeta,
                                               VirusOrder_zeta = VirusOrder_zeta, row_ID = row_ID, arg=arg)}


write.csv(viral_host_sharing_G_df, paste("output/viral_host_sharing_dG_df", array_id, ".csv", sep=""))



