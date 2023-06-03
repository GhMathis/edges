library(tidyverse)
library(doParallel) # for parallel processing
library(parallel) # to detect core

source("R/functions.R")

clover = read.csv("data/clover.csv", stringsAsFactors = T)
trefle = read.csv("data/trefle.csv", stringsAsFactors = T)
##### 1000 interaction for testing /!\to be removed/!\
#set.seed(1000)
#ID = sample(nrow(clover), 1000)
#clover = clover[ID, ]

#trefle = trefle[trefle$host %in% clover$Host & trefle$virus %in% clover$Virus ,]
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
for(h in as.character(levels(trefle$virus))){
  VO = clover$VirusOrder[which(clover$Virus == h)[1]]
  trefle$VirusOrder[trefle$virus == h]= as.character(VO)
}

uni_ntw_trefle = matrix.associations.uni(Virus = trefle$virus, Host = trefle$host)
G_trefle = communicability(uni_ntw_trefle)


trefle$HostOrder = as.character(trefle$HostOrder)
trefle$VirusOrder = as.character(trefle$VirusOrder)

ntw_base = list(uni_ntw_trefle = uni_ntw_trefle,
                G_trefle = G_trefle,
                trefle=trefle)

communi.func<- function(virus,host,HostOrder_zeta, VirusOrder_zeta, row_ID,arg){
  
  print(row_ID)
  
  with(arg,{
    
    ID_virus = which(row.names(uni_ntw_trefle) == virus)  
    ID_host = which(colnames(uni_ntw_trefle) == host) 
    
    uni_ntw_trefle[ID_virus,ID_host] = 0 # change the interaction (perturbation == zeta)
    uni_ntw_trefle[ID_host,ID_virus] = 0 
    
    G_zeta_trefle = communicability(uni_ntw_trefle)
    
    uni_ntw_trefle[ID_virus,ID_host] = 1  # remove the change
    uni_ntw_trefle[ID_host,ID_virus] = 1
    
    
    ### Compute all z score for all group of Host and Virus order.
    # HOST
    # communicability for ALL orders out of the modification and for the order with the modification
    
    zeta_host_order = trefle$host[trefle$HostOrder == HostOrder_zeta]
    zeta_host_extraorder= trefle$host[trefle$HostOrder != HostOrder_zeta]
    
    ID_intraOrder_h =  which(colnames(uni_ntw_trefle) %in% unique(zeta_host_order))

    ID_extraOrder_h =  which(colnames(uni_ntw_trefle) %in% unique(zeta_host_extraorder))
    ID_host = c(ID_intraOrder_h,ID_extraOrder_h)
    #compute G and G zeta means
    g_zeta_intra_order_h = mean(G_zeta_trefle[ID_intraOrder_h, ID_intraOrder_h])
    g_zeta_extra_order_h = mean(G_zeta_trefle[ID_extraOrder_h, ID_extraOrder_h])
    
    g_intra_order_h = mean(G_trefle[ID_intraOrder_h, ID_intraOrder_h])
    g_extra_order_h = mean(G_trefle[ID_extraOrder_h, ID_extraOrder_h])
    
    # mean(G - G zeta) 
    g_d_intra_order_h = mean(G_trefle[ID_intraOrder_h,ID_intraOrder_h]-
      G_zeta_trefle[ID_intraOrder_h, ID_intraOrder_h])
    g_d_extra_order_h = mean(G_trefle[ID_extraOrder_h,ID_extraOrder_h]-
                               G_zeta_trefle[ID_extraOrder_h, ID_extraOrder_h])
    # VIRUS
    VirusOrder_zeta = trefle$VirusOrder[1]
    # communicability for ALL orders out of the modification and for the order with the modification
    zeta_virus_order = trefle$virus[trefle$VirusOrder == VirusOrder_zeta]
    zeta_virus_extraorder= trefle$virus[trefle$VirusOrder != VirusOrder_zeta]
    
    ID_intraOrder_v =  which(colnames(uni_ntw_trefle) %in% unique(zeta_virus_order))
    
    ID_extraOrder_v =  which(colnames(uni_ntw_trefle) %in% unique(zeta_virus_extraorder))
    ID_virus = c(ID_intraOrder_v,ID_extraOrder_v)
    
    #compute G and g zeta means
    g_zeta_intra_order_v = mean(G_zeta_trefle[ID_intraOrder_v, ID_intraOrder_v])
    g_zeta_extra_order_v = mean(G_zeta_trefle[ID_extraOrder_v, ID_extraOrder_v])
    
    g_intra_order_v = mean(G_trefle[ID_intraOrder_v, ID_intraOrder_v])
    g_extra_order_v = mean(G_trefle[ID_extraOrder_v, ID_extraOrder_v])
    
    
    g_d_intra_order_v = mean(G_zeta_trefle[ID_intraOrder_v, ID_intraOrder_v]-
                               G_trefle[ID_intraOrder_v,ID_intraOrder_v])
    g_d_extra_order_v = mean(G_zeta_trefle[ID_extraOrder_v, ID_extraOrder_v] -
                               G_trefle[ID_extraOrder_v,ID_extraOrder_v])
    
    
    df = data.frame(g_zeta_intra_order_h = g_zeta_intra_order_h,
                    g_zeta_extra_order_h = g_zeta_extra_order_h,
                    g_intra_order_h = g_intra_order_h,
                    g_extra_order_h = g_extra_order_h,
                    g_d_intra_order_h = g_d_intra_order_h,
                    g_d_extra_order_h = g_d_extra_order_h,
                    g_zeta_intra_order_v= g_zeta_intra_order_v,
                    g_zeta_extra_order_v = g_zeta_extra_order_v,
                    g_intra_order_v = g_intra_order_v,
                    g_extra_order_v = g_extra_order_v,
                    g_d_intra_order_v = g_d_intra_order_v,
                    g_d_extra_order_v = g_d_extra_order_v,
                    HostOrder_zeta = HostOrder_zeta)
    ##
    
    
    df$mean_dG_h = mean(G_trefle[ID_host, ID_host]-
                         G_zeta_trefle[ID_host, ID_host])
    df$mean_dG_v = mean(G_trefle[ID_virus,ID_virus]-
                          G_zeta_trefle[ID_virus,ID_virus])
    
    df$sd_dG_h = sd(G_trefle[ID_host, ID_host]
                   -G_zeta_trefle[ID_host, ID_host])
    
    df$sd_dG_v = sd(G_zeta_trefle[ID_virus,ID_virus]-
                     G_trefle[ID_virus,ID_virus])
    
    df$X = row_ID
    
    return(df)
    
  })
  
}

print(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
array_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n_array = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))

n_iteration = nrow(trefle)

n = n_iteration %/% n_array
r = n_iteration %% n_array

if(array_id == 1){
  iter = c(1:(n*(array_id)-1))
  print(head(iter,n = 1L))
  print(tail(iter,n = 1L))
}else if(array_id != n_array){
  iter = c((n*array_id):(n*(array_id+1)-1))
  print(head(iter,n = 1L))
  print(tail(iter,n = 1L))
}else if(array_id == n_array){
  iter = c((n*array_id):n_iteration)
  print(head(iter,n = 1L))
  print(tail(iter,n = 1L))
}


ncores =detectCores()
print(ncores)
registerDoParallel(cores=ncores-1)# Shows the number of Parallel Workers to be used
getDoParWorkers()# number of actual workers
viral_host_sharing_G_df = foreach(virus = trefle$virus[iter], host = trefle$host[iter],
                                                                HostOrder_zeta = trefle$HostOrder[iter],
                                                                VirusOrder_zeta = trefle$VirusOrder[iter],
                                                                row_ID =iter,
                                                                arg = lapply(1:length(iter), function (x) ntw_base),
                                                                .combine = rbind, .verbose = T) %dopar%
                                 {communi.func(virus = virus, host = host, HostOrder_zeta = HostOrder_zeta,
                                               VirusOrder_zeta = VirusOrder_zeta, row_ID = row_ID, arg=arg)}


write.csv(viral_host_sharing_G_df, paste("output/viral_host_sharing_dG_df", array_id, ".csv", sep=""))

#host_sharing_G_df = viral_host_sharing_G_df%>%
#  select(g_zeta_intra_order_h,g_zeta_extra_order_h,g_intra_order_h, g_extra_order_h,
#         mean_dG_h,sd_dG_h, HostOrder_zeta,g_d_extra_order_h,g_d_intra_order_h)
#
#
#host_sharing_G_df$score_intra2 =((host_sharing_G_df$g_d_intra_order_h )- #- host_sharing_G_df$g_intra_order_h 
#  (host_sharing_G_df$mean_dG_h ))/#- mean(G_trefle[c(ID_intraOrder_h,ID_extraOrder_h), c(ID_intraOrder_h,ID_extraOrder_h)])
# ( host_sharing_G_df$sd_dG_h )

#host_sharing_G_df$score_extra2 =(( host_sharing_G_df$g_d_extra_order_h )- #- host_sharing_G_df$g_extra_order_h
#  (host_sharing_G_df$mean_dG_h ))/#- mean(G_zeta_trefle[c(ID_intraOrder_h,ID_extraOrder_h), c(ID_intraOrder_h,ID_extraOrder_h)])
# (host_sharing_G_df$sd_dG_h)




