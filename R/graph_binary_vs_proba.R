library(tidyverse)
library(ggpubr)
importance_binary_proba_df = read.csv("output/importance_binary_proba_df.csv", header = T)

main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20),
        legend.title.align=0.5,
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=28))
host_max = importance_binary_proba_df%>%
  count(host)%>%
  arrange(desc(n))%>%
  slice(1:5)%>%
  pull(host)
subset_gg_hostmax = importance_binary_proba_df%>%
  filter(host %in% host_max)
gg_G = ggplot(importance_binary_proba_df)+
  geom_abline(intercept= 0, slope = 1, col = "red", linetype =2, cex =1)+
  geom_point(aes(G_pq_trefle_nrmlz, G_pq_trefle_prob_nrmlz),
             alpha = 0.1,cex =3)+
  geom_point(data=subset_gg_hostmax,aes(G_pq_trefle_nrmlz, G_pq_trefle_prob_nrmlz,
                                   col = host),
             alpha = 0.4,cex =3)+
  main_theme
gg_I = ggplot(importance_binary_proba_df)+
  geom_abline(intercept= 0, slope = 1, col = "red", linetype =2, cex =1)+
  geom_point(aes(importance_trefle_nrmlz, importance_trefle_prob_nrmlz),
             alpha = 0.1,cex =3)+
  geom_point(data=subset_gg_hostmax,aes(importance_trefle_nrmlz, importance_trefle_prob_nrmlz,
                 col =host),
             alpha = 0.4,cex =3)+
  main_theme
ggarrange(gg_G,gg_I)
