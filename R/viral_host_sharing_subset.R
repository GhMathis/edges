library(tidyverse)

temp = read.csv("output/viral_host_sharing_dG_df1.csv", header = T)
for(i in 2:127){
  temp = rbind(temp, read.csv(paste("output/viral_host_sharing_dG_df",i,".csv",sep = ""), header = T))
}
str(temp)
temp = temp[,-1]
write.csv(temp,"output/viral_host_sharing_dG_df.csv")

viral_host_sharing_dG_df = read.csv("output/viral_host_sharing_dG_df.csv", header = T,
                                   stringsAsFactors = T)
trefle = read.csv("data/trefle.csv", header = T,
                                    stringsAsFactors = T)
viral_host_sharing_dG_df = viral_host_sharing_dG_df[,-1]
viral_host_sharing_dG_df$host = trefle$host[1:8959]
viral_host_sharing_dG_df$virus = trefle$virus[1:8959]
str(viral_host_sharing_dG_df)

host_subset = viral_host_sharing_dG_df%>%
  select(ends_with("_h")| c("X","host","HostOrder_zeta"))
str(host_subset)

host_subset$score_intra = ((host_subset$g_d_intra_order_h)- 
                                   (host_subset$mean_dG_h))/(host_subset$sd_dG_h)

host_subset$score_extra = ((host_subset$g_d_extra_order_h )-
                                   (host_subset$mean_dG_h))/ (host_subset$sd_dG_h)
host_subset = host_subset%>%
  pivot_longer(c(score_intra,score_extra), names_to ="score_location", values_to = "score")

write.csv(host_subset, "output/intra_vs_extra_host_subset.csv")
