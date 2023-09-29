
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

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

metadata <- read.table("broilers_metadata.txt", header=TRUE, sep="\t")
str(metadata)

metadata$Day_cat  <- as.factor(metadata$Day_cat)

metadata$Day_cat <- factor(metadata$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

evenness = read_qza("evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("ID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("ID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("ID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("ID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\

simpson = read_qza("simpson_vector.qza")
simpson<-simpson$data %>% rownames_to_column("ID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\

simpson_e = read_qza("simpson_e_vector.qza")
simpson_e<-simpson_e$data %>% rownames_to_column("ID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\

#Make alpha diversity tables
alpha_diversity <- merge(observed_features, shannon, by.x = "ID", by.y = "ID")
alpha_diversity <- merge(alpha_diversity, evenness, by.x = "ID", by.y = "ID")
alpha_diversity <- merge(alpha_diversity, faith_pd, by.x = "ID", by.y = "ID")
alpha_diversity <- merge(alpha_diversity, simpson, by.x = "ID", by.y = "ID")
alpha_diversity <- merge(alpha_diversity, simpson_e, by.x = "ID", by.y = "ID")
meta <- merge(metadata, alpha_diversity, by.x = "ID", by.y = "ID")

str(meta)
meta$observed_features  <- as.numeric(meta$observed_features)

#Save the alpha diversity
write.table(meta,file = "broiler_alpha.txt",quote = F,sep = '\t', row.names = T, col.names = T)


#Mixed model
#Evenness
m1_evenness <- mixed(pielou_evenness ~ Day_cat + Animal + (Day_cat||Animal), data = meta, method = "KR",
                control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_evenness)

anova(m1_evenness)
plot(m1_evenness$full_model)
qqnorm(residuals(m1_evenness$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_evenness<- emmeans(m1_evenness, "Day_cat")
update(pairs(emm_1_evenness), by = NULL, adjust = "holm")
with(meta, shapiro.test(pielou_evenness))

meta$Day_cat <- factor(meta$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

afex_plot(m1_evenness, x = "Day_cat", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Evenness", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/chicken_evenness", ".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#OTU
m1_otu <- mixed(observed_features ~ Day_cat + Animal + (Day_cat||Animal), data = meta, method = "KR",
                control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_otu)

anova(m1_otu)  

plot(m1_otu$full_model)
qqnorm(residuals(m1_otu$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_otu<- emmeans(m1_otu, "Day_cat")
emm_1_otu
update(pairs(emm_1_otu), by = NULL, adjust = "holm")
with(meta, shapiro.test(observed_features))

meta$Day_cat <- factor(meta$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

afex_plot(m1_otu, x = "Day_cat", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Observed ASVs", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/chicken_OTU", ".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#Faith_pd
m1_faith <- mixed(faith_pd ~ Day_cat + Animal + (Day_cat||Animal), data = meta, method = "KR",
                control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_faith)

anova(m1_faith)

plot(m1_faith$full_model)
qqnorm(residuals(m1_faith$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_faith<- emmeans(m1_faith, "Day_cat")
emm_1_faith
update(pairs(emm_1_faith), by = NULL, adjust = "holm")
with(meta, shapiro.test(faith_pd))

meta$Day_cat <- factor(meta$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

afex_plot(m1_faith, x = "Day_cat", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Faith's Phylogenetic Diversity ", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/chicken_faith", ".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

meta <- mutate(meta, faith_log = log10(faith_pd + 1))

#Transformed Faith_pd
T1_faith <- mixed(faith_log ~ Day_cat + Animal + (Day_cat||Animal), data = meta, method = "KR",
                  control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(T1_faith)

anova(T1_faith) 

plot(T1_faith$full_model)
qqnorm(residuals(T1_faith$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_t1_faith<- emmeans(T1_faith, "Day_cat")
emm_t1_faith
update(pairs(emm_t1_faith), by = NULL, adjust = "holm")
with(meta, shapiro.test(faith_log))

meta$Day_cat <- factor(meta$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

afex_plot(T1_faith, x = "Day_cat", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Tansformed Faith's Phylogenetic Diversity ", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/T_chicken_faith", ".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#Shannon
m1_shannon <- mixed(shannon_entropy ~ Day_cat + Animal + (Day_cat||Animal), data = meta, method = "KR",
                  control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_shannon)

anova(m1_shannon) 

plot(m1_shannon$full_model)
qqnorm(residuals(m1_shannon$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_shannon<- emmeans(m1_shannon, "Day_cat")
emm_1_shannon
update(pairs(emm_1_shannon), by = NULL, adjust = "holm")
with(meta, shapiro.test(shannon_entropy))

meta$Day_cat <- factor(meta$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

afex_plot(m1_shannon, x = "Day_cat", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Shannon ", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/chicken_shannon", ".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Simpson
m1_simpson <- mixed(simpson ~ Day_cat + Animal + (Day_cat||Animal), data = meta, method = "KR",
                  control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_simpson)

anova(m1_simpson) 

plot(m1_simpson$full_model)
qqnorm(residuals(m1_simpson$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_simpson<- emmeans(m1_simpson, "Day_cat")
emm_1_simpson
update(pairs(emm_1_simpson), by = NULL, adjust = "holm")
with(meta, shapiro.test(simpson))

meta$Day_cat <- factor(meta$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

afex_plot(m1_simpson, x = "Day_cat", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Simpson ", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/chicken_simpson", ".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#Simpson_evenness
m1_simpson_e <- mixed(simpson_e ~ Day_cat + Animal + (Day_cat||Animal), data = meta, method = "KR",
                    control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_simpson_e)

anova(m1_simpson_e)

plot(m1_simpson_e$full_model)
qqnorm(residuals(m1_simpson_e$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_simpson_e<- emmeans(m1_simpson_e, "Day_cat")
emm_1_simpson_e
update(pairs(emm_1_simpson_e), by = NULL, adjust = "holm")
with(meta, shapiro.test(simpson_e))

meta$Day_cat <- factor(meta$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

afex_plot(m1_simpson_e, x = "Day_cat", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Simpson Evenness", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/chicken_simpson_e", ".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Beta diversity
#Statistics
bc_dist_qza<-read_qza("bray_curtis_distance_matrix.qza")

bc_dist <- bc_dist_qza$data %>% as.matrix() %>%
  as.data.frame(row.names = "SampleID")

meta_dist <- merge (metadata, bc_dist, by=0, all.x = F)
row.names(meta_dist) <- meta_dist$sample.id

meta_dist$Day_cat  <- as.factor(meta_dist$Day_cat)

#bc_dist2 <- bc_dist_qza$data %>% as.matrix()
bc_dist2 <- bc_dist_qza$data
bc_dist2 <- as.matrix(bc_dist2)
#bc_dist2 <- as.dist(bc_dist2)
str(bc_dist2)

library(vegan)
broileer_bc_adnois <- adonis2(as.dist(bc_dist2) ~ meta_dist$Day_cat,
                              strata = meta_dist$Animal)
broileer_bc_adnois2 <- as.matrix(broileer_bc_adnois)
str(broileer_bc_adnois)
broileer_bc_adnois$aov.tab$`Pr(>F)`

bd_broileer_bc <- betadisper(bc_dist2, meta_dist$Day_cat)
anova(bd_broileer_bc)
permutest(bd_broileer_bc)

library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pairwise.adonis2(bc_dist2 ~ Day_cat, data = meta_dist, strata = 'Animal')
pairwise.adonis2(as.dist(bc_dist2) ~ Day_cat, data = meta_dist, strata = 'Animal')


#Plot
bc_PCoA<-read_qza("bray_curtis_pcoa_results.qza")

body_colors <- c("Black", "#fdd49e", "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000")

metadata$Day_cat <- factor(metadata$Day_cat, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

bc_meta <- bc_PCoA$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "ID"))

write.table(bc_meta,file = "BC_broiler.txt",quote = F,sep = '\t', row.names = T, col.names = T)

my_column <- "Day_cat"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Day")
#ggsave(paste0("output/BC-basic_", my_column,"-chicken.pdf"), height=2, width=3, device="pdf") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "Day_cat"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 5) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Day")
ggsave(paste0("output/BC-ellipse_", my_column,"-chicken.pdf"), height=2.5, width=4, device="pdf") # save a PDF 3 inches by 4 inches

#unweighted IniFrac
Uni_PCoA<-read_qza("unweighted_unifrac_pcoa_results.qza")

Uni_meta <- Uni_PCoA$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "ID"))

write.table(Uni_meta,file = "Uni_broiler.txt",quote = F,sep = '\t', row.names = T, col.names = T)

ggplot(Uni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab(paste0("PC1 (", round(100*Uni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Uni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Day")
#ggsave(paste0("output/Un-basic_", my_column,"-chicken.pdf"), height=2, width=3, device="pdf") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Uni_meta,mean)
colnames(centroids)[1] <- "Day_cat"

ggplot(Uni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 5) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Uni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Uni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Day")
ggsave(paste0("output/Uni-ellipse_", my_column,"-chicken.pdf"), height=2.5, width=4, device="pdf") # save a PDF 3 inches by 4 inches

#Weighted UniFrac
Wuni_PCoA<-read_qza("weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "ID"))

write.table(Wuni_meta,file = "Wuni_broiler.txt",quote = F,sep = '\t', row.names = T, col.names = T)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Day")
#ggsave(paste0("output/WU-basic_", my_column,"-chicken.pdf"), height=2, width=3, device="pdf") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)
colnames(centroids)[1] <- "Day_cat"

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 5) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Day")
ggsave(paste0("output/Wuni-ellipse_", my_column,"-chicken.pdf"), height=2.5, width=4, device="pdf") # save a PDF 3 inches by 4 inches


##Qiime2r method of reading in the taxonomy files
taxonomy<-read_qza("taxonomy.qza")
head(taxonomy$data)

tax.clean<-parse_taxonomy(taxonomy$data)
head(tax.clean)


tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}



#################################################################
##Taxa barplot
#################################################################

physeq <- qza_to_phyloseq(
  features="rarefied_table.qza",
  tree="rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "broilers_metadata.txt"
)

metadata <- read.table("broilers_metadata.txt", header=TRUE, sep="\t")
str(metadata)

metadata$Day  <- as.factor(metadata$Day)

metadata$Day  <- as.factor(metadata$Day)
row.names(metadata) <- metadata[,1]

metadata$Day <- factor(metadata$Day, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

#First get the OTU table from physeq
physeq_otu_table <- data.frame(otu_table(physeq), check.names = F)

tax.clean = tax.clean[row.names(tax.clean) %in% rownames(physeq_otu_table),]
metadata.filtered = metadata[row.names(metadata) %in% colnames(physeq_otu_table),]

#Assign as variables to be feed into phyloseq
OTU.physeq = otu_table(as.matrix(physeq_otu_table), taxa_are_rows=TRUE)

#our edited and formatted taxonomy table from the top of this script
tax.physeq = tax_table(as.matrix(tax.clean))    
meta.physeq = sample_data(metadata.filtered)

#We then merge these into an object of class phyloseq.

physeq_bar_plot = phyloseq(OTU.physeq, tax.physeq, meta.physeq)



# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Phylum", "Family", "Genus")
my_column <- "Day"  #this is the metadata column that we will use in the taxa barplot

#rm(taxa.summary)

abund_filter <- 0.05  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$body.site.ord = factor(physeq_meta_filtered$Day, c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("output/", ml, "BarPlot_", my_column, "-chicken.png"), height = 5, width = 5)
}

#Random Forest
# Set working directory
setwd("~/")

getwd()
rm(list = ls()) # clears R memory
package_list <- c("randomForest","ggplot2","pheatmap","vegan","dplyr","magrittr",
                  "scales","grid","reshape2","phyloseq","knitr","qiime2R","tidyr",
                  "naniar","ggpubr","devtools","fansi","pheatmap")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

library(ggplot2)
library(vegan)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(qiime2R)
library(tidyr) #for separate function
library(naniar)# for replace all function
library(ggpubr)
library(pheatmap)

#metadata
metadata <- read.table("broilers_metadata.txt", header=TRUE, sep="\t")
str(metadata)

#metadata$Day  <- as.factor(metadata$Day)
metadata$ID = as.factor(metadata$ID)
rownames(metadata) <- metadata[,1]

#otu_table
otu_table = read.table("feature-rarefied_table-l6.tsv", header = TRUE)


#training
training = metadata[metadata$Animal %in% c("B1","B2","B3","B4", "B5", "B6", "B7", "B8"),] 
#metadata[,1]
#metadata[1,]
idx = rownames(training) %in% colnames(otu_table)
training = training[idx,]

training_otu = otu_table[, rownames(training)]

set.seed(316)
rf = randomForest(t(training_otu), training$Day, importance=TRUE, proximity=TRUE, ntree = 10000)
print(rf)

#features cross validation
set.seed(316) 
#rfcv：Random Forest Cross Validation
result = rfcv(t(training_otu), training$Day, cv.fold=8)
result$error.cv
#plot
p1=with(result, plot(n.var, error.cv, log="x", type="o", lwd=2,xlab="n.var.seed=316"))

cross <- data.frame(result$error.cv)
library(tibble)
cross <- tibble::rownames_to_column(cross, "Number_OTU")
cross = cross[order(1:8, decreasing = T),]
cross = cross[order(cross[,2], decreasing = T),]
head(cross)
write.table(cross,file = "newcross-l6.txt",quote = F,sep = '\t', row.names = T, col.names = T)
cross2 = read.table("newcross-l6.txt", header=T, sep="\t") 


x = ggplot(cross, aes(Number_OTU, result.error.cv, group = 1)) + 
  geom_line()
x =  x + geom_vline(xintercept=30, linetype="dotted", 
                 color = "blue", linewidth=1) +
    labs(y="Cross validation error rate",
         x ="Number of ASVs") +
  theme_classic()
ggsave(paste0("output/Cross_validation_rare-l6", "-pig.png"), x, height = 4, width = 4)
  
imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp)
write.table(imp,file = "importance_genus_rare-16.txt",quote = F,sep = '\t', row.names = T, col.names = T)
#visulisation simple
varImpPlot(rf, main = "Top 30 - Feature OTU importance-16",n.var = 30, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray" )

imp = read.table("importance_genus_rare-16.txt", header=T, row.names= 1, sep="\t") 
imp = head(imp, n=30)
imp = imp[order(1:30,decreasing = T),]

imp$temp = gsub("d__Bacteria;p__","",rownames(imp),perl=TRUE) 

imp$phylum = gsub(";.+","",imp$temp, perl=TRUE)

imp$phylum = gsub("[\\[\\]]+","",imp$phylum,perl=TRUE) 

imp$class = gsub("[\\w\\[\\];_]+;c__","",imp$temp,perl=TRUE)  
imp$class = gsub("[\\[\\]]+","",imp$class,perl=TRUE)

imp$class=factor(imp$class,levels = imp$class)

imp$genus = gsub("[\\w\\[\\];_]+;g__","",imp$temp,perl=TRUE)  
imp$genus = gsub("[\\[\\]]+","",imp$genus,perl=TRUE)

write.table(imp,file = "imp_genus_rare-16.txt",quote = F,sep = '\t', row.names = T, col.names = T)
imp = read.table("imp_genus_rare-16.txt", header=T, row.names= 1, sep="\t") 

imp$genus=factor(imp$genus,levels = imp$genus)
imp = imp[order(1:30,decreasing = T),]

p2=ggplot(data = imp, mapping = aes(x=genus,y=X.IncMSE,fill=phylum)) + 
  geom_bar(stat="identity")+coord_flip()+theme_bw() +
  labs(y = "Increase in mean square error") +
  labs(x = "Genus in increasing importance to accuracy of model") +
  theme(axis.text.x = element_text(color = "black", size = 15), axis.text.y = element_text(color = "black", size = 15)) +
  theme(axis.title.x = element_text(color = "black", size = 14, , face="bold"), axis.title.y = element_text(color = "black", size = 14, face="bold")) 
p2 #importance of taxon in prediction
ggsave(paste0("output/Important_OTUs_RF_rare-l6a", ".pdf"), p2, height=8, width=12, device="pdf") # save a PDF 8 inches by 12 inches


training_abu = training_otu[rownames(imp),]
rownames(training_abu)=imp[rownames(training_abu),"genus"]

p3=pheatmap(training_abu, scale = "row")
p3 #sample vs taxon


####
sampFile = as.data.frame(training$Day,row.names = row.names(training))
colnames(sampFile)[1] = "group"
mat_t = t(training_abu)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
otu_norm_group = do.call(rbind, mat_mean)[-1,]
colnames(otu_norm_group) = mat_mean$group

p4=pheatmap(otu_norm_group,scale="row",cluster_cols = F, cluster_rows = T)
p4 #time vs. taxon
ggsave(paste0("output/Heat_map_rare_Day-l6", "-chicken.png"), p4, height = 8, width = 16)

bak=otu_norm_group
otu_norm_group = otu_norm_group[as.character(imp$genus),] # 按初始排序
for (i in 1:length(rownames(otu_norm_group))) {
  #  i=1
  x=as.data.frame(sort(otu_norm_group[i,],decreasing = T))
  imp[i,"order"]=rownames(x)[1]
}
imp$order2 =  as.numeric(gsub("B1","",imp$order,perl=TRUE) )
imp$order2 =  as.numeric(gsub("B2","",imp$order,perl=TRUE) )
imp$order2 =  as.numeric(gsub("B3","",imp$order,perl=TRUE) )
imp$order2 =  as.numeric(gsub("B4","",imp$order,perl=TRUE) )
imp$order2 =  as.numeric(gsub("B5","",imp$order,perl=TRUE) )
imp$order2 =  as.numeric(gsub("B6","",imp$order,perl=TRUE) )
imp$order2 =  as.numeric(gsub("B7","",imp$order,perl=TRUE) )
imp$order2 =  as.numeric(gsub("B8","",imp$order,perl=TRUE) )
taxonomy = arrange(imp, desc(order2), genus)

otu_norm_group1 = otu_norm_group[match(taxonomy$genus,rownames(otu_norm_group)),] # 按初始排序

p5=pheatmap(otu_norm_group1,scale="row",cluster_cols = F, cluster_rows = F, fontsize_row = 15, fontsize_col = 15)
p5 #time vs taxon 
ggsave(paste0("output/Heat_map_rare_Day_taxon-l6b", "-chicken.png"), p5, height = 8, width = 7)


#using training set to evalue

train.p = predict(rf, type = "response")
df = data.frame(observed = training$Day, predict = train.p)
write.table(df,file = "train_predict-l6.txt",quote = F,sep = '\t', row.names = T, col.names = T)
cor = cor.test(df[,1], df[,2], method = "spearman") # spearman or pearson
m = lm(observed ~ predict, df)
m
summary(m)

p6 = ggplot(df, aes(observed,predict)) +
  geom_point() + geom_jitter()+
  geom_smooth(method = "loess") +
  labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
  theme_bw()
p6
ggsave(paste0("output/Training_prediction_rare-l6", "-chicken.png"), p6, height = 4, width = 4)

test_set = metadata[metadata$Animal%in% c("B9","B10"),]
idx = rownames(test_set) %in% colnames(otu_table)
test_set = test_set[idx,]
test_otu = otu_table[, rownames(test_set)]   
test.p = predict(rf, t(test_otu), type = "response")
df = data.frame(observed = test_set$Day, predict = test.p)
write.table(df,file = "test_predict-l6.txt",quote = F,sep = '\t', row.names = T, col.names = T)
cor = cor.test(df[,1], df[,2], method = "spearman")
m1 = lm(observed ~ predict, df)
m1
summary(m1)
p7 = ggplot(df, aes(observed,predict)) +
  geom_point() + geom_jitter() +
  geom_smooth(method = "loess") +
  labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m1)$r.squared, digits = 3) , sep = "")) +
  theme_bw()
p7
ggsave(paste0("output/Test_prediction_rare-l6", "-chicken.png"), p7, height = 4, width = 4)


#PICRUST
# Set working directory
setwd("~/")

metadata <- read.table("broilers_metadata.txt", header=TRUE, sep="\t")
str(metadata)

metadata$Day <- factor(metadata$Day, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

observed_features = read_qza("observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("ID")

meta <- merge(metadata, observed_features, by.x = "ID", by.y = "ID")

str(meta)
meta$observed_features  <- as.numeric(meta$observed_features)

#OTU
m1_otu <- mixed(observed_features ~ Day + Animal + (Day||Animal), data = meta, method = "KR",
                control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_otu)

anova(m1_otu) ##No significant 
plot(m1_otu$full_model)
qqnorm(residuals(m1_otu$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_otu<- emmeans(m1_otu, "Day")
emm_1_otu
update(pairs(emm_1_otu), by = NULL, adjust = "holm")
with(meta, shapiro.test(observed_features))

afex_plot(m1_otu, x = "Day", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Observed ASVs", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/chicken_picrust_OTU", ".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Beta diversity
bc_PCoA<-read_qza("bray_curtis_pcoa_results.qza")

body_colors <- c("#fff7ec", "#fdd49e", "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000")

metadata$Day <- factor(metadata$Day, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

bc_meta <- bc_PCoA$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "ID"))

write.table(bc_meta,file = "BC_piglet.txt",quote = F,sep = '\t', row.names = T, col.names = T)

my_column <- "Day"

str(bc_meta)

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  #theme_bw() +
  #theme(text=element_text(family = "Garamond"))+
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Day")
#ggsave(paste0("output/BC-basic_", my_column,"-pig.pdf"), height=2, width=3, device="pdf") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "Day"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 5) +
  #theme_classic() +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Day")
ggsave(paste0("output/chicken_picrust_BC-ellipse_", my_column,"-pig.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

# BUGBASE
# Set working directory
setwd("~/bugbase")

metadata <- read.table("broilers_metadata.txt", header=TRUE, sep="\t")
str(metadata)

metadata$Day <- factor(metadata$Day, levels = c("d0", "d1", "d2", "d3", "d7", "d14", "d21"))

phenotypes <- read.table("broilers_predictions.txt", header=TRUE, sep="\t")
str(phenotypes)

meta_phenotypes <- merge(metadata, phenotypes, by.x = "ID", by.y = "ID")
str(meta_phenotypes)

#Anaerobic
m1_Anaerobic <- mixed(Anaerobic ~ Day + Animal + (Day||Animal), data = meta_phenotypes, method = "KR",
                      control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_Anaerobic)

anova(m1_Anaerobic) 
plot(m1_Anaerobic$full_model)
qqnorm(residuals(m1_Anaerobic$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_Anaerobic<- emmeans(m1_Anaerobic, "Day")
emm_1_Anaerobic
update(pairs(emm_1_Anaerobic), by = NULL, adjust = "holm")
with(meta_phenotypes, shapiro.test(Anaerobic))

afex_plot(m1_Anaerobic, x = "Day", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Relative Abundance of Anaerobes", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/broiler_bugbase_Anaerobic", ".pdf"), height=4, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

meta_phenotypes <- mutate(meta_phenotypes, Anaerobic_log = log10(Anaerobic))

# Log Transformed Anaerobic
T1_Anaerobic <- mixed(Anaerobic_log ~ Day + Animal + (Day||Animal), data = meta_phenotypes, method = "KR",
                      control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(T1_Anaerobic)

anova(T1_Anaerobic) 
plot(T1_Anaerobic$full_model)
qqnorm(residuals(T1_Anaerobic$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_t1_Anaerobic<- emmeans(T1_Anaerobic, "Day")
emm_t1_Anaerobic
update(pairs(emm_t1_Anaerobic), by = NULL, adjust = "holm")
with(meta_phenotypes, shapiro.test(Anaerobic_log))

afex_plot(T1_Anaerobic, x = "Day", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Relative Abundance of Log Anaerobes", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/T_broiler_bugbase_Anaerobic", ".pdf"), height=4, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#Aerobic
m1_Aerobic <- mixed(Aerobic ~ Day + Animal + (Day||Animal), data = meta_phenotypes, method = "KR",
                    control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_Aerobic)

anova(m1_Aerobic) 
plot(m1_Aerobic$full_model)
qqnorm(residuals(m1_Aerobic$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_Aerobic<- emmeans(m1_Aerobic, "Day")
emm_1_Aerobic
update(pairs(emm_1_Aerobic), by = NULL, adjust = "holm")
with(meta_phenotypes, shapiro.test(Aerobic))

afex_plot(m1_Aerobic, x = "Day", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Relative Abundance of Aerobes", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/broiler_bugbase_Aerobic", ".pdf"), height=4, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#Facultatively Anaerobic
m1_F_Anaerobic <- mixed(Facultatively_Anaerobic ~ Day + Animal + (Day||Animal), data = meta_phenotypes, method = "KR",
                        control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_F_Anaerobic)

anova(m1_F_Anaerobic) 
plot(m1_F_Anaerobic$full_model)
qqnorm(residuals(m1_F_Anaerobic$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_F_Anaerobic<- emmeans(m1_F_Anaerobic, "Day")
emm_1_F_Anaerobic
update(pairs(emm_1_F_Anaerobic), by = NULL, adjust = "holm")
with(meta_phenotypes, shapiro.test(Facultatively_Anaerobic))

afex_plot(m1_F_Anaerobic, x = "Day", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Relative Abundance of Facultative Anaerobes", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/broiler_bugbase_F_Anaerobic", ".pdf"), height=4, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Forms_Biofilms
m1_Biofilms <- mixed(Forms_Biofilms ~ Day + Animal + (Day||Animal), data = meta_phenotypes, method = "KR",
                     control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_Biofilms)

anova(m1_Biofilms) 
plot(m1_Biofilms$full_model)
qqnorm(residuals(m1_Biofilms$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_Biofilms<- emmeans(m1_Biofilms, "Day")
emm_1_Biofilms
update(pairs(emm_1_Biofilms), by = NULL, adjust = "holm")
with(meta_phenotypes, shapiro.test(Forms_Biofilms_log))

afex_plot(m1_Biofilms, x = "Day", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Relative Abundance of Biofilm Formers", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/broiler_bugbase_Biofilms", ".pdf"), height=4, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Gram_Negative
m1_Gram_Negative <- mixed(Gram_Negative ~ Day + Animal + (Day||Animal), data = meta_phenotypes, method = "KR",
                          control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_Gram_Negative)

anova(m1_Gram_Negative) 
plot(m1_Gram_Negative$full_model)
qqnorm(residuals(m1_Gram_Negative$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_Gram_Negative<- emmeans(m1_Gram_Negative, "Day")
emm_1_Gram_Negative
update(pairs(emm_1_Gram_Negative), by = NULL, adjust = "holm")
with(meta_phenotypes, shapiro.test(Gram_Negative))

afex_plot(m1_Gram_Negative, x = "Day", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Relative Abundance of Gram_Negative", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/broiler_bugbase_Gram_Negative", ".pdf"), height=4, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Gram_Positive
m1_Gram_Positive <- mixed(Gram_Positive ~ Day + Animal + (Day||Animal), data = meta_phenotypes, method = "KR",
                          control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_Gram_Positive)

anova(m1_Gram_Positive) 
plot(m1_Gram_Positive$full_model)
qqnorm(residuals(m1_Gram_Positive$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_Gram_Positive<- emmeans(m1_Gram_Positive, "Day")
emm_1_Gram_Positive
update(pairs(emm_1_Gram_Positive), by = NULL, adjust = "holm")
with(meta_phenotypes, shapiro.test(Gram_Positive))

afex_plot(m1_Gram_Positive, x = "Day", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Relative Abundance of Gram Positive", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/broiler_bugbase_Gram_Positive", ".pdf"), height=4, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Contains_Mobile_Elements
m1_Contains_Mobile_Elements <- mixed(Contains_Mobile_Elements ~ Day + Animal + (Day||Animal), data = meta_phenotypes, method = "KR",
                                     control = lmerControl(optCtrl = list(maxfun=1e6)), expand_re=TRUE)

summary(m1_Contains_Mobile_Elements)

anova(m1_Contains_Mobile_Elements) 
plot(m1_Contains_Mobile_Elements$full_model)
qqnorm(residuals(m1_Contains_Mobile_Elements$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1_Contains_Mobile_Elements<- emmeans(m1_Contains_Mobile_Elements, "Day")
emm_1_Contains_Mobile_Elements
update(pairs(emm_1_Contains_Mobile_Elements), by = NULL, adjust = "holm")
with(meta_phenotypes, shapiro.test(Contains_Mobile_Elements))

afex_plot(m1_Contains_Mobile_Elements, x = "Day", id = "ID", error = "between", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") + 
  labs(y = "Relative Abundance of Bacteria with Mobile Elements", x = "Day") +
  theme(legend.position="none") 

ggsave(paste0("output/broiler_bugbase_Contains_Mobile_Elements", ".pdf"), height=4, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches
