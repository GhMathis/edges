library(tidyverse)
library(doParallel) # for parallel processing
library(parallel) # to detect core

source("R/functions.R")

communicability.lenght.path<- function(A, deepness){
  spectra = eigen(A)
  d = spectra$values
  P = spectra$vectors
  out = 0
  k=0
  
  G_list = list()
  k=1
  for (k in deepness){
    #d = D**k
    if(k == "inf"){
      
      out = communicability(A, spectra =spectra )
    }else{
      out = out + 1/factorial(k)*(P %*%  diag(d**k) %*% solve(P)) # A + A²
      # A + P%*%D
     
    }
    G_list = append(G_list, list(out) )
  }
  
  return(G_list)
}
#test = communicability.lenght.path(matrix(c(1,1,1,1), ncol = 2, nrow = 2), c(2,4,6,8))

clover = read.csv("data/clover.csv", stringsAsFactors = T)
trefle = read.csv("data/trefle.csv", stringsAsFactors = T)
deepness = list(2,6,10,20,50,"inf")
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
G_trefle = communicability.lenght.path(uni_ntw_trefle, deepness = deepness)

# trefle$IDmatrix_v = sapply(trefle$virus, function(x) which(row.names(uni_ntw_trefle) == x)) 
# trefle$IDmatrix_h = sapply(trefle$host, function(x) which(row.names(uni_ntw_trefle) == x)) 
trefle$HostOrder = as.character(trefle$HostOrder)
trefle$VirusOrder = as.character(trefle$VirusOrder)
HO = trefle%>%
  count(HostOrder)%>%
  arrange(desc(n))%>%pull(HostOrder)
VO = trefle%>%
  count(VirusOrder)%>%
  arrange(desc(n))%>%pull(VirusOrder)
trefle_main_HostOrder = trefle %>%
  filter(HostOrder %in% HO[1:15])
trefle_main_VirusOrder = trefle %>%
  filter(VirusOrder %in% VO[1:15])

ntw_base = list(uni_ntw_trefle = uni_ntw_trefle,
                G_trefle = G_trefle,
                trefle_main_HostOrder = trefle_main_HostOrder,
                trefle_main_VirusOrder = trefle_main_VirusOrder,
                deepness = deepness)

communi.func<- function(virus,host,HostOrder_zeta, VirusOrder_zeta, row_ID,arg){
  
  print(row_ID)
  
  with(arg,{
    ID_virus = which(row.names(uni_ntw_trefle) == virus)  
    ID_host = which(colnames(uni_ntw_trefle) == host) 
    
    ###preturbed adjacency matrix
    uni_ntw_trefle[ID_virus,ID_host] = 0 # change the interaction (perturbation == zeta)
    uni_ntw_trefle[ID_host,ID_virus] = 0 
    
    # preturbed communicability matrix
    G_zeta_trefle = communicability.lenght.path(uni_ntw_trefle, deepness)
    
    # adjacency matrix
    uni_ntw_trefle[ID_virus,ID_host] = 1  # remove the change
    uni_ntw_trefle[ID_host,ID_virus] = 1
    
    # df initialization
    df = data.frame(HostOrder_zeta = character(),
                    VirusOrder_zeta = character(),
                    HostOrder = character(),
               host = character(),
               score_host_group = numeric(),
               VirusOrder = character(),
               virus = character(),
               score_virus_group = numeric(),
               deepness = integer(), X = integer())
    # df_h = data.frame(HostOrder_zeta = character(),
    #                 host = character(),
    #                 score_host_group = numeric(),
    #                 deepness = integer(), X = integer())
    # 
    # df_v =data.frame( VirusOrder_zeta = character(),
    #                 virus = character(),
    #                 score_virus_group = numeric(),
    #                 deepness = integer(), X = integer())
    i = 6
    for(d in deepness){
      i = i+1
      ##### Compute all z score for all groups of Host and Virus order.
      
      ### HOST     
      
      # global effect of the perturbation 
      
      host_host_ID =  which(colnames(uni_ntw_trefle) %in% trefle$host)

      mean_dG_h = mean(G_trefle[[i]][host_host_ID, host_host_ID]-
                         G_zeta_trefle[[i]][host_host_ID, host_host_ID])
      
      sd_dG_h = sd(G_trefle[[i]][host_host_ID, host_host_ID]
                   -G_zeta_trefle[[i]][host_host_ID, host_host_ID])
      score_host_group  = c()
      for(HO in unique(trefle_main_HostOrder$HostOrder)){
        
        #
        host_group = trefle$host[trefle$HostOrder == HO]
        ID_host_group =  which(colnames(uni_ntw_trefle) %in% host_group)
      
        # mean(G - G zeta) 
        g_d_host_group = mean(G_trefle[[i]][ID_host_group,ID_host_group]-
                                   G_zeta_trefle[[i]][ID_host_group, ID_host_group])
        
        score_host_group = c(score_host_group, ((g_d_host_group)- 
                                                  (mean_dG_h))/(sd_dG_h))
      }
      HostOrder = unique(trefle_main_HostOrder$HostOrder)
      
      ### VIRUS
      # global effect of the perturbation 
      
      virus_virus_ID =  which(colnames(uni_ntw_trefle) %in% trefle$virus)
      
      mean_dG_h = mean(G_trefle[[i]][virus_virus_ID, virus_virus_ID]-
                         G_zeta_trefle[[i]][virus_virus_ID, virus_virus_ID])
      sd_dG_h = sd(G_trefle[[i]][virus_virus_ID, virus_virus_ID]
                   -G_zeta_trefle[[i]][virus_virus_ID, virus_virus_ID])
      score_virus_group  =c()
      for(VO in unique(trefle_main_VirusOrder$VirusOrder)){
        
        #
        virus_group = trefle$virus[trefle$VirusOrder == VO]
        ID_virus_group =  which(colnames(uni_ntw_trefle) %in% virus_group)
        
        # mean(G - G zeta) 
        g_d_virus_group = mean(G_trefle[[i]][ID_virus_group,ID_virus_group]-
                                G_zeta_trefle[[i]][ID_virus_group, ID_virus_group])
        
        score_virus_group = c(score_virus_group, ((g_d_virus_group)- 
                                                  (mean_dG_h))/(sd_dG_h))
      }
      VirusOrder = unique(trefle_main_HostOrder$HostOrder)
      # df_h =rbind(df_h,data.frame(HostOrder_zeta = HostOrder_zeta,
      #                   host = host,
      #                   score_host_group = score_host_group,
      #                   deepness = d, X = row_ID))
      # 
      df =rbind(df,data.frame(HostOrder_zeta = HostOrder_zeta,
                              VirusOrder_zeta = VirusOrder_zeta,
                              HostOrder = HostOrder,
                        host = host,
                        score_host_group = score_host_group,
                        VirusOrder = VirusOrder,
                        virus = virus,
                        score_virus_group = score_virus_group,
                        deepness = d, X = row_ID))
      
    }

    return(df)
    
  })
  
}
### subset of interaction
set.seed(1000)


trefle$X = 1:nrow(trefle)
iter = trefle%>%
  group_by(HostOrder)%>%
  slice_sample(n=1000)%>%
  filter(HostOrder %in% unique(trefle_main_HostOrder$HostOrder))%>%
  pull(X)
str(as.factor(trefle[iter,]$VirusOrder))
#iter = sample(nrow(trefle), 5000)
n_iteration = length(iter)
###

print(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
array_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n_array = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT")) 


n = n_iteration %/% n_array
r = n_iteration %% n_array

if(array_id == 1){
 ID = c(1:(n*(array_id)-1))
 print(head(ID,n = 1L))
 print(tail(ID,n = 1L))
}else{
  ID = c((n*(array_id-1)):(n*(array_id)-1))
 print(head(ID,n = 1L))
 print(tail(ID,n = 1L))
}
if(array_id == n_array){
  ID = c(ID,c((n*array_id):n_iteration))
 print(head(ID,n = 1L))
 print(tail(ID,n = 1L))
}
iter = iter[ID]
ncores =detectCores()
print(ncores)
registerDoParallel(cores=ncores-1)# Shows the number of Parallel Workers to be used
getDoParWorkers()# number of actual workers


viral_host_sharing_extended_df = foreach(virus = trefle$virus[iter], host = trefle$host[iter],
                                  HostOrder_zeta = trefle$HostOrder,
                                  VirusOrder_zeta = trefle$VirusOrder,
                                  row_ID = as.integer(row.names(trefle)[iter]),
                                  arg = lapply(1:length(iter), function (x) ntw_base),
                                  .combine = rbind,.verbose = T) %dopar%
  {communi.func(virus = virus, host = host, HostOrder_zeta = HostOrder_zeta,
                VirusOrder_zeta = VirusOrder_zeta, row_ID = row_ID, arg=arg)}

#write.csv(viral_host_sharing_extended_df, "output/viral_host_sharing_extended_df.csv")
write.csv(viral_host_sharing_extended_df, paste("output/viral_host_sharing_extended",
                                          "_df", array_id, ".csv", sep=""))

