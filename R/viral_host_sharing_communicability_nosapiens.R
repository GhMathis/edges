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


trefle = trefle %>%
  filter(host != "Homo sapiens")

trefle$HostOrder= "_"
trefle$VirusOrder= "_"
trefle$host = as.factor(trefle$host)
trefle$virus = as.factor(trefle$virus)
# attribue host order
for(h in as.character(levels(trefle$hos))){
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

trefle$IDmatrix_v = sapply(trefle$virus, function(x) which(row.names(uni_ntw_trefle) == x)) 
trefle$IDmatrix_h = sapply(trefle$host, function(x) which(row.names(uni_ntw_trefle) == x)) 
trefle$HostOrder = as.character(trefle$HostOrder)
trefle$VirusOrder = as.character(trefle$VirusOrder)

trefle_main_HostOrder = trefle %>%
  filter(HostOrder %in% c("Rodentia", "Primates", "Cetartiodactyla", "Carnivora",
                          "Chiroptera", "Perissodactyla", "Lagomorpha", "Diprotodontia",
                           "Didelphimorphia", "Eulipotyphla"))
trefle_main_VirusOrder = trefle %>%
  filter(VirusOrder %in% c("Bunyavirales", "Mononegavirales", "Herpesvirales", "Amarillovirales",
                            "Picornavirales", "Reovirales", "Martellivirales","Chitovirales",
                            "Piccovirales", "Zurhausenvirales"))

ntw_base = list(uni_ntw_trefle = uni_ntw_trefle,
                G_trefle = G_trefle,
                trefle_main_HostOrder = trefle_main_HostOrder,
                trefle_main_VirusOrder = trefle_main_VirusOrder)
virus = trefle$virus[1]
host = trefle$host[1]
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
    HO ="Didelphimorphia"
    
    ### Compute all z score for all group of Host and Virus order.
   
    score_intra_h =c()
    score_intra_v =c()
    score_global_h =c()
    score_global_v =c()
    # HOST
    
    # unique score for ALL orders out of the modification
    zeta_host_order= trefle_main_HostOrder$host[trefle_main_HostOrder$HostOrder == HostOrder_zeta]
    ID_zetaOrder_h =  which(colnames(uni_ntw_trefle) %in% zeta_host_order)
    
    #compute z-score for all host order exept host order with the zeta modification 
    score_extra_order_h =(mean(G_zeta_trefle[-ID_zetaOrder_h, -ID_zetaOrder_h]) -
                            mean(G_trefle[-ID_zetaOrder_h, -ID_zetaOrder_h]))/sd(G_trefle)#[-ID_zetaOrder_h, -ID_zetaOrder_h])
    # unique score for EACH orders 
    for(HO in unique(trefle_main_HostOrder$HostOrder)){
      all_host_per_order = trefle_main_HostOrder$host[trefle_main_HostOrder$HostOrder == HO]
     
      # removed link host order
      ID_HostOrder_intra = which(colnames(uni_ntw_trefle) %in% all_host_per_order) 
    
      ### group mean modification compaire to INTRA (itself)
      temp = (mean(G_zeta_trefle[ID_HostOrder_intra,ID_HostOrder_intra]) - mean(G_trefle[ID_HostOrder_intra,ID_HostOrder_intra]))/sd(G_trefle)#[ID_HostOrder_intra,ID_HostOrder_intra])

      score_intra_h = c(score_intra_h,temp)
      ### group mean modification compaire to the global communicability
      temp = (mean(G_zeta_trefle[ID_HostOrder_intra,ID_HostOrder_intra]) - mean(G_zeta_trefle))/sd(G_zeta_trefle)
      score_global_h = c(score_global_h,temp)
    }
     all_host_per_order = trefle_main_HostOrder$host[trefle_main_HostOrder$HostOrder == HO]
    # VIRUS
    
   # unique score for ALL orders out of the modification
    zeta_virus_order= trefle_main_VirusOrder$virus[trefle_main_VirusOrder$VirusOrder == VirusOrder_zeta]
    ID_zetaOrder_v =  which(colnames(uni_ntw_trefle) %in% zeta_virus_order)
    
    #compute z-score for all host order exept host order with the zeta modification 
    score_extra_order_v =(mean(G_zeta_trefle[-ID_zetaOrder_v, -ID_zetaOrder_v]) -
                            mean(G_trefle[-ID_zetaOrder_v, -ID_zetaOrder_v]))/sd(G_trefle)#[-ID_zetaOrder_v, -ID_zetaOrder_v])
    # unique score for EACH orders 
    for(VO in unique(trefle_main_VirusOrder$VirusOrder)){
      all_virus_per_order = trefle_main_VirusOrder$virus[trefle_main_VirusOrder$VirusOrder == VO]

      
      # removed link virus order
      ID_VirusOrder_intra = which(colnames(uni_ntw_trefle) %in% all_virus_per_order) 
    
      ## group mean modification compaire to INTRA (itself)
      temp = (mean(G_zeta_trefle[ID_VirusOrder_intra,ID_VirusOrder_intra]) - mean(G_trefle[ID_VirusOrder_intra,ID_VirusOrder_intra]))/sd(G_trefle)#[ID_VirusOrder_intra,ID_VirusOrder_intra])
      
      score_intra_v = c(score_intra_v,temp)
      
      ## group mean modification compaire to the global communicability
      temp = (mean(G_zeta_trefle[ID_VirusOrder_intra,ID_VirusOrder_intra]) - mean(G_zeta_trefle))/sd(G_zeta_trefle)
      score_global_v = c(score_global_v,temp)
    }
    df = data.frame(HostOrder = unique(trefle_main_HostOrder$HostOrder),
                    score_intra_h = score_intra_h, host_zeta = host,
                    VirusOrder = unique(trefle_main_VirusOrder$VirusOrder),
                    score_intra_v = score_intra_v, virus_zeta = virus,
                    score_extra_order_h = score_extra_order_h,
                    score_extra_order_v = score_extra_order_v,
                    HostOrder_zeta = HostOrder_zeta,
                    VirusOrder_zeta = VirusOrder_zeta,
                    score_global_h = score_global_h,
                    score_global_v = score_global_v)
    ##
  

    df$complet_score =(mean(G_zeta_trefle)- mean(G_trefle))/sd(G_trefle)
    
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
registerDoParallel(cores=(ncores-1))# Shows the number of Parallel Workers to be used
getDoParWorkers()# number of actual workers


viral_host_sharing_G_df = foreach(virus = trefle$virus[iter], host = trefle$host[iter],
                                   HostOrder_zeta = trefle$HostOrder,
                                   VirusOrder_zeta = trefle$VirusOrder,
                                   row_ID = row.names(trefle)[iter],
                                   arg = lapply(1:length(iter), function (x) ntw_base),
                                   .combine = rbind, .verbose = T) %dopar%
  {communi.func(virus = virus, host = host, HostOrder_zeta = HostOrder_zeta,
                VirusOrder_zeta = VirusOrder_zeta, row_ID = row_ID, arg=arg)}

write.csv(viral_host_sharing_G_df, paste("output/viral_host_sharing_G_df", array_id, ".csv", sep=""))


