library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(car)
library(ggforce)
main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=20),
        axis.text.y = element_text(colour = "black", size=20),
        legend.title = element_text(colour = "black", size=22),
        legend.title.align=0,
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=22))

###### Viral sharing #######
host_subset = read.csv("output/intra_vs_extra_host_subset.csv", header = T,
                                   stringsAsFactors = T)

host_subset$score_location = as.factor(host_subset$score_location)
levels(host_subset$score_location) = c("Extra Ordre", "Intra Ordre")
str(host_subset)
pal = brewer.pal(n =9, name = "OrRd")

levels(host_subset$HostOrder_zeta)
main_order = host_subset%>%
  count(HostOrder_zeta)%>%
  arrange(desc(n))%>%
  slice(1:10)%>%
  pull(HostOrder_zeta)
recap2 = host_subset%>%
  filter(HostOrder_zeta %in% unique(main_order))%>%
  group_by(score_location,HostOrder_zeta)%>%
  summarise(N    =length(score),
            mean_score = mean(score, na.rm=TRUE),
            sd_score = sd(score, na.rm=TRUE),
            se_score = sd_score / sqrt(N))
host_subset%>%
  filter(HostOrder_zeta %in% unique(main_order))%>%
ggplot()+
  geom_sina(aes(HostOrder_zeta, score, col = score_location),scale = "width", position=position_dodge(0.5),
            maxwidth = 0.4, cex =2, alpha = 0.05)+
  geom_errorbar(data = recap2, aes(x= HostOrder_zeta, group = score_location,
                                  ymin=mean_score-sd_score, ymax=mean_score+sd_score),
                cex =0.8, width=.3, position = position_dodge(0.5)) +
  geom_point(data = recap2, aes(x= HostOrder_zeta, y=mean_score,
                               group = score_location), cex =3,
             position=position_dodge(0.5)) +
  color_palette(pal[c(8,6,4)])+
  labs(col =c("Calcule du partage viral"), x = "Ordre", y = "Partage Viral")+
  main_theme+
  guides(colour = guide_legend(override.aes = list(alpha = 1, cex =3)))+
  theme(axis.text.x = element_text(angle =45, hjust = 1),
        #egend.justification=c(1,0),
        legend.position = "top")

recap = intra_vs_extra_sharing_nosapiens%>%
  filter(score_location != "Intra + Extra Ordre")%>%
  group_by(score_location,HostOrder_zeta)%>%
  summarise(N    =length(z_score),
            mean_zscore = mean(z_score, na.rm=TRUE),
            sd_zscore = sd(z_score, na.rm=TRUE),
            se_zscore = sd_zscore / sqrt(N))

intra_vs_extra_sharing_nosapiens%>%
  filter(score_location != "Intra + Extra Ordre")%>%
  ggplot()+
  geom_sina(aes(HostOrder_zeta, z_score, col = score_location),scale = "width", position=position_dodge(0.5),
            maxwidth = 0.4, cex =2, alpha = 0.05)+
  geom_errorbar(data = recap, aes(x= HostOrder_zeta, group = score_location,
                                  ymin=mean_zscore-sd_zscore, ymax=mean_zscore+sd_zscore),
                cex =0.8, width=.3, position = position_dodge(0.5)) +
  geom_point(data = recap, aes(x= HostOrder_zeta, y=mean_zscore,
                               group = score_location), cex =3,
              position=position_dodge(0.5)) +
  color_palette(pal[c(8,6,4)])+
  labs(col =c("Calcule du partage viral"), x = "Ordre", y = "Partage Viral")+
  main_theme+
  guides(colour = guide_legend(override.aes = list(alpha = 1, cex =3)))+
  theme(axis.text.x = element_text(angle =45, hjust = 1),
        #egend.justification=c(1,0),
        legend.position = "top")


mod_v_sharing = lm(z_score~HostOrder_zeta*score_location, data = intra_vs_extra_sharing_nosapiens%>%
                     filter(score_location != "Intra + Extra Ordre"))

par(mfrow = c(2,2))
plot(mod_v_sharing)
hist(mod_v_sharing$residuals)
summary(mod_v_sharing)
Anova(mod_v_sharing, type = 2)
TukeyHSD(aov(mod_v_sharing))
###### Trefle vs Clover #######

## importance and communicability score for observed associassion (trefle and clover)
importance_df = read.csv("output/importance_df.csv")
summary(importance_df)
length(importance_df$status)
G_nrmlz_subset = importance_df %>%
  select(c(virus, host, X, G_pq_clover_nrmlz, G_pq_trefle_nrmlz))%>%
  pivot_longer(-c(virus, host, X),values_to ="G_nrmlz")
I_nrmlz_subset = importance_df %>%
  select(c(virus, host, X, importance_trefle_nrmlz, importance_clover_nrmlz))%>%
  pivot_longer(-c(virus, host, X),values_to ="I_nrmlz")

#### data whitout sapiens
## importance and communicability score for all associassion (trefle)
importance_trefle_nosapiens_df = read.csv("output/importance_trefle_nosapiens_df.csv")
importance_nosapiens_df = read.csv("output/importance_nosapiens_df.csv")


hist_G = ggplot(importance_df)+
  geom_density(aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz))+
  #geom_smooth(aes(G_pq_clover_nrmlz,G_pq_trefle_nrmlz), method = "lm", color = "black", linetype = 2, cex = 2)+
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  labs(fill ="Compte du nombre \n associations", x = "Communicabilité CLOVER", y= "Communicabilité TREFLE")+
  main_theme+
  theme(legend.title.align = 0.5,
        legend.direction = "vertical",
        legend.box.just = "center")+
  coord_fixed()


hist_G = ggplot(G_nrmlz_subset)+
  geom_density(aes(G_nrmlz, fill = name),alpha = 0.3, stat = "density")+
  fill_palette("pastel1")+
  main_theme+
  scale_fill_discrete(name = "Réseaux", labels = c("CLOVER", "TREFLE"))+
  labs(x ="Communicabilité (G)", y= "Densité")

hist_I = ggplot(I_nrmlz_subset)+
  geom_density(aes(I_nrmlz, fill = name),alpha = 0.3, stat = "density")+
  fill_palette("pastel1")+
  scale_fill_discrete(name = "Réseaux", labels = c("CLOVER", "TREFLE"))+
  labs(x ="Importance (I)", y= "Densité") +
  main_theme

subset_gg_hostmax = importance_df%>%
  filter(host %in% c(importance_df%>%
                       count(host)%>%
                       arrange(desc(n))%>%
                       slice(1:10)%>%pull(host)))
subset_gg_hostrest = importance_df%>%
  filter(!(host %in% c(importance_df%>%
                       count(host)%>%
                       arrange(desc(n))%>%
                       slice(1:10)%>%pull(host))))
subset_gg_G2 = subset(importance_df, G_pq_trefle_nrmlz > 0.75 |
                        G_pq_clover_nrmlz > 0.25)
gg_G2  = ggplot(subset_gg_hostrest)+
  #geom_bin_2d(aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz))+
  geom_abline(intercept= 0, slope = 1, col = "red", linetype =2, cex =1)+
  geom_point(aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz), alpha = 0.1,cex =2)+
  geom_point(data=subset_gg_hostmax,
             aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz,col = host), alpha = 0.4, cex =2)+
  
  color_palette("set3") +
  
  # geom_text_repel(data=subset_gg_G2,
  #                 aes(G_pq_clover_nrmlz,G_pq_trefle_nrmlz, label = paste(host, "(", virus, ")", sep = "")),
  #                 size=4,
  #                 box.padding = unit(0.5, "lines"),
  #                 max.overlaps =20)+
  labs(x ="Communicabilité (G colver)", y= "Communicabilité (G trefle)", col = "Hôtes")+
  main_theme

subset_gg_I2 = subset(importance_df, importance_trefle_nrmlz > 0.009 |
                        (importance_trefle_nrmlz > 0.008 & importance_clover_nrmlz > 0.03) |
                        importance_clover_nrmlz > 0.043)

gg_I2 = ggplot(subset_gg_hostrest)+
  geom_abline(intercept= 0, slope = 1, col = "red", linetype =2, cex =1)+
  geom_point(aes(importance_clover_nrmlz, importance_trefle_nrmlz), alpha = 0.1,cex =2)+
  geom_point(data=subset_gg_hostmax,
             aes(importance_clover_nrmlz, importance_trefle_nrmlz,col = host), alpha = 0.4, cex =2)+
  color_palette("set3") +
  # geom_text_repel(data=subset_gg_I2,
  #                 aes(importance_clover_nrmlz,importance_trefle_nrmlz, label = paste(host, "(", virus, ")", sep = "")),
  #                 size=4,
  #                 box.padding = unit(0.5, "lines"),
  #                 max.overlaps = 20)+
  # 
  labs(x ="Importance (I colver)", y= "Importance (I trefle)", col = "Hôtes")+
  main_theme
gg_trefle = ggplot(importance_df)+
  geom_point(aes(G_pq_trefle_nrmlz, importance_trefle_nrmlz))+
  geom_abline(intercept= 0, slope = 1, col = "blue")+
  labs(x ="G trefle", y= "I trefle")+
  main_theme
ggarrange(hist_G,hist_I,gg_G2,gg_I2, ncol =2,nrow = 2)



###### Importance vs chance of association (trefle only) #######
imputed_association = read.csv("data/imputed_associations.csv")
#### data with sapiens
## importance and communicability score for imputed associassion (trefle only)
importance_trefle_df = read.csv("output/importance_unshared_df.csv")
importance_trefle_df$status = "imputed associations"

## importance and communicability score for observed associassion (trefle and clover)
importance_df = read.csv("output/importance_df.csv")
importance_df$status = "observed associations"

#### data whitout sapiens
# data whitout sapiens

importance_trefle_nosapiens_df = read.csv("output/importance_trefle_nosapiens_df.csv")
importance_nosapiens_df = read.csv("output/importance_nosapiens_df.csv")

v = importance_trefle_nosapiens_df$virus
h = importance_trefle_nosapiens_df$host

ID = which(paste(v,h) %in% paste(imputed_association$virus, imputed_association$host))
importance_trefle_nosapiens_df$status[ID] = "imputed associations"
importance_trefle_nosapiens_df$status[-ID] = "observed associations"
importance_nosapiens_df[ID,]


# joint to have association proba and importance in same the df
str(imputed_association)
imputed_association = imputed_association%>%
  full_join(importance_trefle_df,by =c("host", "virus"))
str(imputed_association)

# Imputed + observed associations
importance_trefle_df = rbind(importance_trefle_df, importance_df[,names(importance_df) %in% names(importance_trefle_df)])

str(importance_trefle_df)

ggplot(imputed_association)+
  geom_bin_2d(aes(G_pq_trefle_nrmlz,evidence), bins = 100)+
  geom_smooth(aes(G_pq_trefle_nrmlz,evidence), method = "lm", color = "black", linetype = 2, cex = 2)+
  scale_y_log10()+
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  labs(fill ="Compte du nombre \n associations", x = "Communicabilité", y= "evidence (échelle log)")+
  main_theme+
  theme(legend.title.align = 0.5,
        legend.direction = "vertical",
        legend.box.just = "center")

######
# Imputed + observed associations

str(importance_trefle_df)
# imputed + observed and ordering to have coresspndig col names
importance_trefle_df = rbind(importance_trefle_df,
                             importance_df[names(importance_df) %in% names(importance_trefle_df)])
ggplot(importance_trefle_df)+
  #geom_bin_2d(aes(G_pq_clover_nrmlz, G_pq_trefle_nrmlz))+
  geom_point(aes(G_pq_trefle_nrmlz, importance_trefle_nrmlz, col= status), alpha = 0.1,cex =2)+
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  # geom_text_repel(data=subset_gg_G2,
  #                 aes(G_pq_clover_nrmlz,G_pq_trefle_nrmlz, label = paste(host, "(", virus, ")", sep = "")),
  #                 size=4,
  #                 box.padding = unit(0.5, "lines"),
  #                 max.overlaps =20)+
  coord_cartesian(xlim = c(0,0.77))+
  labs(x ="G trefle", y= "I trefle")+
  main_theme
