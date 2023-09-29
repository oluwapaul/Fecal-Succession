library(stringi)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(emmeans)
library(lme4)
library(lattice)
library("latticeExtra")
library(dplyr)
library(summarytools)
library(ggfortify)
library(tidyverse)
library(multcomp)
library(ggbeeswarm)
library(phyloseq)
library(DESeq2)
library(qiime2R)

# Set working directory
setwd("~/")

R.version()

#Combine alpha diversity
combine_alpha = read.table("combine_alpha.txt", header=T, sep="\t") 
combine_alpha$Day_cat <- factor(combine_alpha$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

#combine_alpha <- mutate(type = factor(type, levels = unique(type)))
ggplot(combine_alpha, aes(x=`Day`, y=observed_features, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Observed ASV") +
  theme_classic() +
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("Observed_ASV_by_time2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_alpha, aes(x=`Day`, y=faith_pd, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Faith's Phylogenetic Diversity") +
  theme_classic() +
 # theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("Faith_pd_by_time2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_alpha, aes(x=`Day`, y=pielou_evenness, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Evenness") +
  theme_classic() +
  # theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("Evenness_by_time2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_alpha, aes(x=`Day_cat`, y=shannon_entropy, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Shannon") +
  theme_classic() +
  # theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("Shannon_by_time2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches


#Combine beta diversity (PC1)
combine_beta = read.table("combine_beta.txt", header=T, sep="\t") 
combine_beta$Day_cat <- factor(combine_beta$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

ggplot(combine_beta, aes(x=`Day`, y=BC, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Bray-Curtis PC1") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("BC_by_time.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_beta, aes(x=`Day_cat`, y=BC, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Bray-Curtis PC1") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("BC_by_time2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(combine_beta, aes(x=`Day`, y=UN, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Unweighted UniFrac  PC1") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("UN_by_time.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_beta, aes(x=`Day_cat`, y=UN, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Unweighted UniFrac PC1") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("UN_by_time2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_beta, aes(x=`Day`, y=WU, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Weighted UniFrac  PC1") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("WU_by_time.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_beta, aes(x=`Day_cat`, y=WU, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Weighted UniFrac PC1") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("WU_by_time2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches


#Combine beta diversity2 (Distance from 0)
combine_beta2 = read.table("combine_beta_distance.txt", header=T, sep="\t") 
combine_beta2$Day_cat <- factor(combine_beta2$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

ggplot(combine_beta2, aes(x=`Day`, y=BC, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Bray-Curtis Distance from Day 0") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("BC_by_time_dis.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_beta2, aes(x=`Day_cat`, y=BC, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Bray-Curtis Distance from Day 0") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("BC_by_time_dis2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(combine_beta2, aes(x=`Day`, y=UN, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Unweighted UniFrac  Distance Day 0") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("UN_by_time_dis.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_beta2, aes(x=`Day_cat`, y=UN, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Unweighted UniFrac Distance 0") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("UN_by_time_dis2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_beta2, aes(x=`Day`, y=WU, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Weighted UniFrac  Distance from 0") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("WU_by_time_dis.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(combine_beta2, aes(x=`Day_cat`, y=WU, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Weighted UniFrac Distance from 0") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("WU_by_time_dis2.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches


#Normalized Combine rarefied phylum abundance (Top 5 phylum), numalization was done by dividing by the number of sequences used for rarefaction
norm_combine_phylum = read.table("norm_combine_phylum_abund2.txt", header=T, sep="\t") 
#norm_combine_phylum$Day_cat <- factor(combine_phylum$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))


#Bacteroidetes
ggplot(norm_combine_phylum, aes(x=`Day`, y=Bacteroidetes, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Percentage of Bacteroidetes") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("norm_rare_Bacteroidetes.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

#Firmicutes
ggplot(norm_combine_phylum, aes(x=`Day`, y=Firmicutes, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Percentage of Firmicutes") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("norm_rare_Firmicutes.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

#Proteobacteria
ggplot(norm_combine_phylum, aes(x=`Day`, y=Proteobacteria, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Percentage of Proteobacteria") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("norm_rare_Proteobacteria.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

#Actinobacteria
ggplot(norm_combine_phylum, aes(x=`Day`, y=Actinobacteria, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Percentage of Actinobacteria") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("norm_rare_Actinobacteria.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

#Spirochaetes
ggplot(norm_combine_phylum, aes(x=`Day`, y=Spirochaetes, color=`type`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Percentage of Spirochaetes") +
  theme(text=element_text(family="Garamond", size=14))+
  theme_classic()+
  #theme_q2r() +  # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="type") # use different color scale which is color blind friendly
ggsave("norm_rare_Spirochaetes.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches
