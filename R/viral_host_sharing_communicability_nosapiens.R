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
                          "Chiroptera", "Perissodactyla", "Lagomorpha"))
trefle_main_VirusOrder = trefle %>%
  filter(VirusOrder %in% c("Bunyavirales", "Mononegavirales", "Herpesvirales", "Amarillovirales",
                            "Picornavirales", "Reovirales", "Martellivirales") )

ntw_base = list(uni_ntw_trefle = uni_ntw_trefle,
                G_trefle = G_trefle,
                trefle_main_HostOrder = trefle_main_HostOrder,
                trefle_main_VirusOrder = trefle_main_VirusOrder)

communi.func<- function(virus,host,row_ID,arg){
  
  print(row_ID)
  
  with(arg,{
    
    ID_virus = which(row.names(uni_ntw_trefle) == virus)  
    ID_host = which(colnames(uni_ntw_trefle) == host) 
    uni_ntw_trefle[ID_virus,ID_host] = 0 # change the interaction (perturbation == zeta)
    uni_ntw_trefle[ID_host,ID_virus] = 0 

    G_zeta_trefle = communicability(uni_ntw_trefle)
    
    uni_ntw_trefle[ID_virus,ID_host] = 1  # remove the change
    uni_ntw_trefle[ID_host,ID_virus] = 1
    HO =unique(trefle_main_HostOrder$HostOrder)[1]
    
    ### Compute all z score for all group of Host and Virus order.

    score_intra_h =c()
    score_extra_intra_h = c()
    score_intra_v =c()
    score_extra_intra_v = c()
    
    for(HO in unique(trefle_main_HostOrder$HostOrder)){
      all_host_per_order = trefle_main_HostOrder$host[trefle_main_HostOrder$HostOrder == HO]
      all_host_per_order_extra = trefle_main_HostOrder$host[!(trefle_main_HostOrder$HostOrder == HO)]
      # removed link host order
      ID_HostOrder_intra = which(colnames(uni_ntw_trefle) %in% all_host_per_order) 
      #
      ID_HostOrder_extra = which(colnames(uni_ntw_trefle) %in% all_host_per_order_extra) 
      # ID_HostOrder = ID_HostOrder[ID_HostOrder != ID_host]
      
      ### group mean modification compaire to INTRA (itself)
      temp = (mean(G_zeta_trefle[ID_HostOrder_intra,ID_HostOrder_intra]) - mean(G_trefle[ID_HostOrder_intra,ID_HostOrder_intra]))/sd(G_trefle[ID_HostOrder_intra,ID_HostOrder_intra])

      score_intra_h = c(score_intra_h,temp)
      
      ### EXTRA group modification compare to INTRA groupe modification
      temp = (mean(G_zeta_trefle[ID_HostOrder_extra,ID_HostOrder_extra]) - mean(G_trefle[ID_HostOrder_intra,ID_HostOrder_intra]))/sd(G_trefle[ID_HostOrder_intra,ID_HostOrder_intra])

      score_extra_intra_h = c(score_extra_intra_h,temp)
    }
    temp_v = c()
    VO = unique(trefle_main_VirusOrder$VirusOrder)[2]
    for(VO in unique(trefle_main_VirusOrder$VirusOrder)){
      all_virus_per_order = trefle_main_VirusOrder$virus[trefle_main_VirusOrder$VirusOrder == VO]
      all_virus_per_order_extra = trefle_main_VirusOrder$virus[!(trefle_main_VirusOrder$VirusOrder == VO)]
      
      # removed link virus order
      ID_VirusOrder_intra = which(colnames(uni_ntw_trefle) %in% all_virus_per_order) 
      #
      ID_VirusOrder_extra = which(colnames(uni_ntw_trefle) %in% all_virus_per_order_extra) 
      # ID_VirusOrder = ID_VirusOrder[ID_VirusOrder != ID_virus]
      
      ## group mean modification compaire to INTRA (itself)
      temp = (mean(G_zeta_trefle[ID_VirusOrder_intra,ID_VirusOrder_intra]) - mean(G_trefle[ID_VirusOrder_intra,ID_VirusOrder_intra]))/sd(G_trefle[ID_VirusOrder_intra,ID_VirusOrder_intra])
      
      score_intra_v = c(score_intra_v,temp)
      
      ### EXTRA group modification compare to INTRA groupe modification
      temp = (mean(G_zeta_trefle[ID_VirusOrder_extra,ID_VirusOrder_extra]) - mean(G_trefle[ID_VirusOrder_intra,ID_VirusOrder_intra]))/sd(G_trefle[ID_VirusOrder_intra,ID_VirusOrder_intra])

      score_extra_intra_v = c(score_extra_intra_v,temp)
    }
    df = data.frame(HostOrder = unique(trefle_main_HostOrder$HostOrder),
                    score_extra_intra_h = score_extra_intra_h, 
                    score_intra_h = score_intra_h, host_zeta = host,
                    VirusOrder = unique(trefle_main_VirusOrder$VirusOrder),
                    score_intra_v = score_intra_v, virus_zeta = virus,
                    score_extra_intra_v = score_extra_intra_v)
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


intraOrder_importance_df = foreach(virus = trefle$virus[iter], host = trefle$host[iter],
                             row_ID = row.names(trefle)[iter],
                             arg = lapply(1:length(iter), function (x) ntw_base),
                             .combine = rbind, .verbose = T) %dopar%
  {communi.func(virus = virus, host = host, row_ID = row_ID, arg=arg)}

write.csv(intraOrder_importance_df, paste("output/intraOrder_importance_df", array_id, ".csv", sep=""))


