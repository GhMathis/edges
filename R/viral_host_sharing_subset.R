library(tidyverse)

viral_host_sharing_G_df = read.csv("output/viral_host_sharing2_G_df.csv", header = T,
                                   stringsAsFactors = T)
viral_host_sharing_G_df = viral_host_sharing_G_df[,-1]

viral_host_sharing_G_df$zeta_HostOrder_location = ifelse(as.character(viral_host_sharing_G_df$HostOrder)==
                                                           as.character(viral_host_sharing_G_df$HostOrder_zeta),
                                                         "intra", "extra")
viral_host_sharing_G_df$zeta_VirusOrder_location = ifelse(as.character(viral_host_sharing_G_df$VirusOrder)==
                                                            as.character(viral_host_sharing_G_df$VirusOrder_zeta),
                                                          "intra", "extra")

intra_vs_extra_sharing_subset = viral_host_sharing_G_df%>%
  select(c("HostOrder_zeta","score_intra_h","host_zeta" ,"score_extra_order_h",
           "complet_score", "zeta_HostOrder_location", "X"))%>%
  filter(zeta_HostOrder_location =="intra")%>%
  select(-c( "zeta_HostOrder_location"))%>%
  pivot_longer(c(score_intra_h,score_extra_order_h,complet_score), names_to = "score_location", values_to = "score")
str(intra_vs_extra_sharing_subset)
ID_extra =which(intra_vs_extra_sharing_subset$score_location == "score_extra_order_h")
ID_intra =which(intra_vs_extra_sharing_subset$score_location == "score_intra_h")
ID_global =which(intra_vs_extra_sharing_subset$score_location == "complet_score")
intra_vs_extra_sharing_subset$z_score = NA
intra_vs_extra_sharing_subset$z_score[ID_extra] = (intra_vs_extra_sharing_subset$score[ID_extra]-
                                           mean(intra_vs_extra_sharing_subset$score[ID_global]))/
  sd(intra_vs_extra_sharing_subset$score[ID_global])
intra_vs_extra_sharing_subset$score[ID_extra][1]
intra_vs_extra_sharing_subset$z_score[ID_intra] = (intra_vs_extra_sharing_subset$score[ID_intra]-
                                                     mean(intra_vs_extra_sharing_subset$score[ID_global]))/
  sd(intra_vs_extra_sharing_subset$score[ID_global])
intra_vs_extra_sharing_subset$z_score[ID_global] = (intra_vs_extra_sharing_subset$score[ID_global]-
                                                     mean(intra_vs_extra_sharing_subset$score[ID_global]))/
  sd(intra_vs_extra_sharing_subset$score[ID_global])

write.csv(intra_vs_extra_sharing_subset, "output/intra_vs_extra_sharing_nosapiens.csv")
