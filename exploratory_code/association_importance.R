
library(tidyverse)
library(svd)

setwd("~/Fac/Master_Rennes/stage1/trefle")
prediction = read.csv("hpc/outputs/predictions.csv", stringsAsFactors = T)
str(prediction)
head(prediction)
plot(1/factorial(1:10), (1/2)**(1:10))
clover = read.csv("data/clover.csv", stringsAsFactors = T)
str(clover)
head(clover)
names(clover)

imputed_associations = read.csv("artifacts/imputed_associations.csv", stringsAsFactors = T)
str(imputed_associations)
head(imputed_associations)
geom_
trefle = read.csv("artifacts/trefle.csv", stringsAsFactors = T)
str(trefle)
svd_test =svd(A)
eigen(A)
svd_test$u[1,1]*svd_test$d[1]*t(svd_test$v)[1,1]
D=diag(svd_test$d)
D = cbind(D,0)
A
round(svd_test$u)%*%t(round(svd_test$v))%*%round(svd_test$u)%*%t(round(svd_test$v))
round(svd_test$u %*% diag(svd_test$d)%*%t(svd_test$v))

##### Setup the associations matrix
matrix.associations = function(Virus, Host){
  long_df = data.frame(Virus = Virus, Host = Host) 
  n_virus = length(levels(Virus))
  n_host = length(levels(Host))
  ntw = matrix(0, nrow = n_virus, ncol = n_host)
  
  colnames(ntw) = levels(Host)
  rownames(ntw) = levels(Virus)
  for(v in unique(Virus)){
    temp = long_df %>% filter(Virus == v)
    for(h in unique(temp$Host)){
      ntw[which(rownames(ntw) == v),
                 which(colnames(ntw) == h)] = 1
      
    }
  }
  return(ntw)
}
ntw_clover = matrix.associations(clover$Virus, clover$Host)
sum(ntw_clover)/(ncol(ntw_clover)*nrow(ntw_clover))

ntw_trefle = matrix.associations(trefle$virus, trefle$host)
sum(ntw_trefle)/(ncol(ntw_trefle)*nrow(ntw_trefle))

table(ntw_clover != ntw_trefle)
nrow(imputed_associations)

# Get the matrice of all new imputed associations
test = imputed_associations[subset_imp_asso,]
test = test%>%
  filter(dist_L ==0)
new_imputed_associations <- imputed_associations[!(paste(imputed_associations$virus, imputed_associations$host) %in%
                                      paste(clover$Virus, clover$Host)),]
test[!(paste(test$virus, test$host) %in%
                         paste(trefle$virus, trefle$host)),]


##### compute SVD and random grph dot product
svd.rpd <-function(ntw){

  SVD = propack.svd(ntw, neig= 12)
  trunc_d =  diag(sqrt(SVD$d))
  #diag(trunc_d)[13:length(diag(trunc_d))] = 0

  ### truncate at 12 dim
  L =(SVD$u %*% trunc_d)#[,1:12]
  egein_V = SVD$d
 
  R = trunc_d %*% t(SVD$v)
  t_R = t(R)#[,1:12]

  return(list(L = L,t_R = t_R, egein_V = egein_V))
}

##### subset of the dataset to build the code
# TO REMOVE FOR REAL ANALYSIS
###
# ntw_clover = ntw_clover[1:50,1:50] 
# 
# imputed_associations = imputed_associations %>%
#   filter(host %in% colnames(ntw_clover), virus %in% rownames(ntw_clover))
# str(ntw_clover)

###
# Time test
# 
T1=Sys.time()
for(i in 1:10){
  I0 = svd.rpd(ntw_trefle)

}
T2 = Sys.time()

(T2-T1)*nrow(trefle)/(10*3600)

# This matrix is used as reference for distance
I0 = svd.rpd(ntw_trefle)

base_L = norm(I0$L)
base_R = norm(I0$t_R)
dist_L = vector(length = nrow(imputed_associations))
dist_R = vector(length = nrow(imputed_associations))
dist_elementwise_L = vector(length = nrow(imputed_associations))
dist_elementwise_R = vector(length = nrow(imputed_associations))
subset_imp_asso = sample(1:nrow(imputed_associations), 10000)
str(subset_imp_asso)
subset_imp_asso = sort.int(subset_imp_asso)
svd_data = vector("list", length(subset_imp_asso))
length(svd_data)



for(n in subset_imp_asso){# 1:nrow(imputed_associations)
 
  ID = c(which(rownames(ntw_trefle) == imputed_associations$virus[n]),
               which(colnames(ntw_trefle) == imputed_associations$host[n]))
 
  ntw_trefle[ID[1],ID[2]] = 0 # change the interaction
  I_temp = svd.rpd(ntw_trefle)
  ntw_trefle[ID[1],ID[2]] = 1 # remove the change
  svd_data[[which(subset_imp_asso == n)]] = I_temp # save svd.rpd

  dist_L[n] = base_L - norm(I_temp$L)
  dist_R[n] = base_R - norm(I_temp$t_R)
  
  dist_elementwise_L[n]= sum(colSums(abs(I0$L - I_temp$L)))
  dist_elementwise_R[n]= sum(colSums(abs(I0$t_R - I_temp$t_R)))
}


imputed_associations$dist_L = dist_L
imputed_associations$dist_R = dist_R
imputed_associations$dist_elementwise_L = dist_elementwise_L
imputed_associations$dist_elementwise_R = dist_elementwise_R
str(imputed_associations)
ggplot(imputed_associations[1:58140,])+
  geom_point(aes(dist_L,P))+
  geom_smooth(aes(dist_L,P))
ggplot(imputed_associations[1:58140,])+
  geom_point(aes(dist_R,P))+
  geom_smooth(aes(dist_R,P))
ggplot(imputed_associations[1:58140,])+
  geom_point(aes(dist_elementwise_L,P))+
  geom_smooth(aes(dist_elementwise_L,P))
 
save(svd_data,file = "Mathis/svd_10000.RData")
str(imputed_associations)

imputed_associations[subset_imp_asso,] %>%
  filter(dist_L ==0)

imputed_associations$subset_asso = 0
imputed_associations[subset_imp_asso,]$subset_asso = subset_imp_asso
str(imputed_associations)
write.csv(imputed_associations[subset_imp_asso,], "Mathis/test_10000_asso.csv")
# Other approch egienvector centrality
link.importance<- function(A, deepness){
  
  # scaling value given the first eigen value
  lambda1 = eigen(A)$values[1]
  alpha = 1/lambda1 
  
  #
  D = diag(eigen(A)$values)
  P = eigen(A)$vectors
  d = D
  out = A
  for (deep in 1:deepness){
    d =(alpha**deep)*d*D
    out = out + P %*% d %*% solve(P)
  }
  
  return(out)
}
test =  matrix(c(0,1,
                 1,1), nrow =2, ncol =2)
test = matrix(c(0,0,1,0,
                0,0,1,1,
                1,1,0,0,
                0,1,0,0), nrow =4, ncol =4)
test%*%test

test4=link.importance(test,5)
rowSums(test4)
eigen(test4)

matrix.associations.uni = function(Virus, Host){
  long_df = data.frame(Virus = Virus, Host = Host) 
  n_virus = length(unique(Virus))
  n_host = length(unique(Host))
  ntw = matrix(0, nrow = n_virus+n_host, ncol = n_host+n_virus)
  
  colnames(ntw) = c(unique(Host),unique(Virus))
  rownames(ntw) = c(unique(Host),unique(Virus))
  for(v in unique(Virus)){
    temp = long_df %>% filter(Virus == v)
    for(h in unique(temp$Host)){
      ntw[which(rownames(ntw) == v),
          which(colnames(ntw) == h)] = 1
      ntw[which(rownames(ntw) == h),
          which(colnames(ntw) == v)] = 1
      
    }
  }
  return(ntw)
}
uni_ntw_trefle = matrix.associations.uni(trefle$virus, trefle$host)
str(uni_ntw_trefle)
importance_matrix = link.importance(uni_ntw_trefle,deepness= 5)


imputed_associations$importance = 0
for(n in 1:nrow(imputed_associations)){
  ID = c(which(rownames(importance_matrix) == imputed_associations$virus[n]),
         which(colnames(importance_matrix) == imputed_associations$host[n]))
  imputed_associations$importance[n] = importance_matrix[ID[1],ID[2]]
}
str(imputed_associations)
ggplot(imputed_associations[1:58140,], aes(importance, dist_L))+
  geom_point()
ggplot(imputed_associations[1:58140,], aes(importance, dist_R))+
  geom_point()
ggplot(imputed_associations[1:58140,], aes(n_association, importance))+
  geom_point()
ggplot(imputed_associations[1:58140,], aes(importance,P ))+
  geom_point()+
  geom_smooth()
# Other approch with markov chain

##### Test on test on the importance of a species
# Change all the imputed association of on specie and compaire it's SVD
# to the originale matrix's SVD
# guess that sp with more iteraction will ahve more impact
# but 2nd or 3th order interaction might be also relevante
# test_matrice = matrix(c(0,0,1,1,
#                         1,0,0,1,
#                         0,0,0,1,
#                         1,1,1,0,
#                         0,1,0,0), ncol = 4, nrow =5)
# row.names(test_matrice)= c("a","b","c","d","e")
# colnames(test_matrice) = c( "Macropus robustus", "Cervus nippon",
#                              "Meles meles", "Sylvilagus brasiliensis")
dist_L_row = vector(length = nrow(ntw_trefle))
dist_R_row = vector(length = nrow(ntw_trefle))
n_imputed_association = vector(length = nrow(ntw_trefle))
n_imputed_interaction = c()
T_all = c()

for(v in rownames(ntw_trefle)){
  host_subset =  imputed_associations%>%
    filter(v == imputed_associations$virus)%>%
    pull(host)
  ID_virus = which(row.names(ntw_trefle) == v)  
  ID_host = which(colnames(ntw_trefle) %in% host_subset) 
  
  ntw_trefle[ID_virus,ID_host] = 0 # change the interaction
  I_temp = svd.rpd(ntw_trefle)
  ntw_trefle[ID_virus,ID_host] = 1 # remove the change
  dist_L_row[ID_virus] = base_L - norm(I_temp$L)
  dist_R_row[ID_virus] =base_R - norm(I_temp$t_R)
  
  #imputed association only
  n_imputed_association[ID_virus] = length(host_subset)
}

# all association
(T_all[2] - T_all[1])*1000
n_association = rowSums(ntw_trefle)
n_imputed_association[1:10]
str(n_association[1:10])
par(mfrow = c(1,2))
plot((dist_L_row/n_imputed_association)~n_imputed_association)
plot((dist_L_row/n_imputed_association)~n_association, col = "blue", type="p")

virus_distances = data.frame(n_imputed_association = n_imputed_association,
                             n_association = n_association,
                             dist_L_row = dist_L_row,
                             dist_R_row = dist_R_row)

virus_distances  = virus_distances%>%
  filter( n_imputed_association !=0 )

ggplot(virus_distances,aes(n_imputed_association,
                           dist_L_row/n_imputed_association))+
  geom_point()+
  geom_smooth()
ggplot(virus_distances,aes(n_imputed_association,
                           dist_R_row))+
  geom_point()+
  geom_smooth()
ggplot(virus_distances,aes(n_association,
                           dist_R_row))+
  geom_point()+
  geom_smooth()

##### distance 
str(n_association[names(n_association) %in%
                    imputed_associations$virus])
imputed_associations$n_association = 0
for (v in names(n_association)){
  ID = which(imputed_associations$virus == v)  
  if(length(ID)!=0){
    imputed_associations[ID,]$n_association = n_association[names(n_association) == v]
  }
}
str(imputed_associations)
ggplot(imputed_associations[1:58140,], aes(n_association,dist_R))+
  geom_point()+
  geom_smooth()+
  xlab("nbr of imputed interaction for virus")+
  ylab("distance of R matrix after one interaction modification")
ggplot(imputed_associations[1:58140,], aes(n_association,dist_L))+
  geom_point()+
  geom_smooth()+
  xlab("nbr of imputed interaction for virus")+
  ylab("distance of L matrix after one interaction modification")
str(I0$L)
str(I0$t_R)


