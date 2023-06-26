library(tidyverse)

temp = read.csv("output/viral_host_sharing_dG_df1.csv", header = T)
for(i in 2:127){
  temp = rbind(temp, read.csv(paste("output/viral_host_sharing_dG_df",i,".csv",sep = ""), header = T))
}
temp = temp[,-(names(temp) == "X.1")]

write.csv(temp,"output/viral_host_sharing_dG_df.csv")
