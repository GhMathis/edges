str(viral_host_sharing_extended_df)
host_recap = viral_host_sharing_extended_df%>%
  select(c("HostOrder_zeta","HostOrder","host", "score_host_group","deepness","X"))%>%
  group_by(HostOrder,HostOrder_zeta, deepness)%>%
  summarise(mean_score = mean(score_host_group, na.rm=TRUE),
            sd_score = sd(score_host_group, na.rm=TRUE))

ggplot(host_recap)+
  geom_point(aes(HostOrder,HostOrder_zeta , col = mean_score),size = 5)+
  facet_wrap(~deepness)


temp = viral_host_sharing_extended_df%>%
  select(c("HostOrder_zeta","HostOrder","host", "score_host_group","deepness","X"))%>%
  filter(HostOrder_zeta ==HostOrder )
  group_by(HostOrder,HostOrder_zeta, deepness)%>%
  summarise(mean_score = mean(score_host_group, na.rm=TRUE),
            sd_score = sd(score_host_group, na.rm=TRUE))