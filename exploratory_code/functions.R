library(svd)
normalized = function(x) (x-min(x))/(max(x)-min(x))

##### From row ntw to matrice ntw (bipartie => non_square)
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
##### From row ntw to matrice ntw (unipartie => square)
# give a symetrical square matrice Host+Virus x Host+Virus association
# all Host-Host or Virus-Virus associations are 0 
matrix.associations.uni = function(Virus, Host){
  long_df = data.frame(Virus = Virus, Host = Host) 
  n_virus = length(unique(Virus))
  n_host = length(unique(Host))
  ntw = matrix(0, nrow = n_virus+n_host, ncol = n_host+n_virus) # put host and virus
  # in both colums and rows
  
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
  norm
  return(ntw)
}
##### svd and RDPG function 
svd.rpd <-function(ntw){
  ### truncate at 12 dim
  SVD = propack.svd(ntw, neig= 12)
  trunc_d =  diag(sqrt(SVD$d))
  
  
  
  L =(SVD$u %*% trunc_d)
  egein_V = SVD$d
  
  R = trunc_d %*% t(SVD$v)
  t_R = t(R)#[,1:12]
  
  return(list(L = L,t_R = t_R, egein_V = egein_V))
}
##### Communicability combute with spectra of the A matrix
communicability<- function(A, spectra = NULL){
  if(is.null(spectra)){
    spectra = eigen(A)
  }
  d = spectra$values
  P = spectra$vectors
  G=0
  for(i in 1:length(d)){
    phi = P[,i]/norm(P[,i],"2")
    temp = phi%*%t(phi)*exp(d[i])
    G = G + temp
  }
  rownames(G)= rownames(A)
  colnames(G)= colnames(A)
  return(G)
}
ntw = matrix.associations(trefle$virus, trefle$host)
str(ntw)
communicability_svd <- function(A){
  SVD = propack.svd(ntw)
  U = SVD$u
  V = SVD$v
  str(V)
  cosh(SVD$d[1])
  sinh(SVD$d[1])
  sinh(SVD$d)[1:20]
  d_left =  diag(cosh(SVD$d))
  d_right =  diag(sinh(SVD$d))
  
  G_left = U%*%d_left%*%t(U)
  G_right = U%*%d_right%*%t(V)
  G_svd = G_left - G_right
  str(G_svd)
  rownames(G)= rownames(A)
  colnames(G)= colnames(A)
  return(G)
}