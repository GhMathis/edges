viral_host_sharing_extended_df = read.csv("output/viral_host_sharing_extended_df.csv", header = T)
str(viral_host_sharing_extended_df)
host_recap = viral_host_sharing_extended_df%>%
  select(c("HostOrder_zeta","HostOrder","host", "score_host_group","deepness","X"))%>%
  group_by(HostOrder,HostOrder_zeta, deepness)%>%
  summarise(mean_score = mean(score_host_group, na.rm=TRUE),
            sd_score = sd(score_host_group, na.rm=TRUE))
viral_host_sharing_extended_df$score_host_group
ggplot(host_recap)+
  geom_point(aes(HostOrder_zeta,HostOrder , col = mean_score),size = 10,shape = 15)+
  facet_wrap(~deepness)

max(host_recap$mean_score)
temp = viral_host_sharing_extended_df%>%
  select(c("HostOrder_zeta","HostOrder","host", "score_host_group","deepness","X"))%>%
  filter(HostOrder_zeta ==HostOrder )
  group_by(HostOrder,HostOrder_zeta, deepness)%>%
  summarise(mean_score = mean(score_host_group, na.rm=TRUE),
            sd_score = sd(score_host_group, na.rm=TRUE))