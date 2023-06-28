library(tidyverse)

temp = read.csv("output/importance_binary_proba_df1.csv", header = T)
for(i in 2:100){
  temp = rbind(temp, read.csv(paste("output/importance_binary_proba_df",i,".csv",sep = ""), header = T))
}
importance = read.csv("output/importance_df.csv", header = T)
importance_unshared = read.csv("output/importance_unshared_df.csv", header = T)
str(importance_unshared)
importance_trefle = rbind(importance_unshared%>%
                            select(c("X","virus","host","G_pq_trefle_nrmlz",
                                     "importance_trefle_nrmlz")),
                          importance%>%
                            select(c("X","virus","host","G_pq_trefle_nrmlz",
                                     "importance_trefle_nrmlz")))
str(importance_trefle)
str(temp)
importance_binary_proba_df = temp %>%
  select(-c("virus","host"))%>%
  full_join(importance_trefle, by = "X")
str(importance_binary_proba_df)
write.csv(importance_binary_proba_df,"output/importance_binary_proba_df.csv")
