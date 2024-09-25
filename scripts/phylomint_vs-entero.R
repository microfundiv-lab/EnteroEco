# load libraries
library(data.table)
library(tidyverse)
library(ggsignif)
library(ggrastr)
library(ggpubr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/gem")
metadata = read.csv("../overlap/combined_coloniz-abund_filt.csv")
rownames(metadata) = metadata$Species_rep
metadata = metadata[which(metadata$Family != "Enterobacteriaceae"),]

# parse phylomint results
phylomint.res = read.table("phylomint_gapless.tsv", fill=TRUE, header=TRUE)
phylomint.res$identifier = apply(phylomint.res[, c("A", "B")], 1, function(x) paste(sort(x), collapse="_"))
phylomint.res = phylomint.res[!duplicated(phylomint.res$identifier), ]
phylomint.res$identifier = NULL

# parse results for Enterobacteriaceae
entero.species = read.delim("../metadata/species_entero.txt", header = F, stringsAsFactors = F)
phy.entero = phylomint.res[which((phylomint.res$A %in% metadata$Species_rep & phylomint.res$B %in% entero.species$V1) |
                                   (phylomint.res$B %in% metadata$Species_rep & phylomint.res$A %in% entero.species$V1)),]
phy.entero$Species_rep = ifelse(phy.entero$A %in% entero.species$V1, phy.entero$B, phy.entero$A)
phy.entero$Entero = ifelse(phy.entero$A == phy.entero$Species_rep, phy.entero$B, phy.entero$A)
phy.entero$Distance = 1-(as.numeric(phy.entero$Competition) - as.numeric(phy.entero$Complementarity))
phy.entero = merge(phy.entero, metadata[c("Species_rep", "Classification", "Order", "Family", "Genus")], by = "Species_rep")

# analyse clusters
phy.entero$Cluster = NA
phy.entero$Prod = phy.entero$Competition*phy.entero$Complementarity
phy.entero$Cluster[which(phy.entero$Prod < 0.06)] = "Cluster2"
phy.entero$Cluster = ifelse(is.na(phy.entero$Cluster), "Cluster1", "Cluster2")

# cluster dfs
clust1.df = phy.entero[which(phy.entero$Cluster == "Cluster1"),]
clust2.df = phy.entero[which(phy.entero$Cluster == "Cluster2"),]

# check correlations
clust1.corr = cor.test(clust1.df$Competition, clust1.df$Complementarity)
clust2.corr = cor.test(clust2.df$Competition, clust2.df$Complementarity)

# plot and save
scatter.entero = ggplot(phy.entero, aes(x=Competition, y=Complementarity, colour=Classification)) +
  geom_point_rast(size=0.8) +
  scale_colour_manual(values = c("Co-excluder" = 'steelblue', "Co-colonizer" = "tomato")) +
  geom_smooth(aes(colour=Cluster), method="lm", se=TRUE) +
  ylab("Metabolic complementarity") +
  xlab("Metabolic competition") +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  theme_classic() +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12))

plot.split = ggplot(phy.entero, aes(x=Classification, y= Distance, fill = Classification)) +
    geom_point_rast(colour = "darkgrey", size=0.5, position = position_jitterdodge(jitter.width = 0.5)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    facet_wrap(~ Cluster) +
    theme_classic() +
    scale_x_discrete(limits=c("Co-excluder", "Co-colonizer")) +
    scale_fill_manual(values = c("Co-excluder" = 'steelblue', "Co-colonizer" = "tomato")) +
    geom_signif(
      comparisons = list(c("Co-excluder", "Co-colonizer")),
      map_signif_level = FALSE) +
    ylab("Metabolic distance") +
    theme(strip.background = element_blank()) +
    theme(strip.text = element_text(size=14)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size=14)) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.text.x = element_text(size=12))

# combine cluster plots
ggarrange(scatter.entero, plot.split, common.legend=TRUE, labels=c("a", "b"), font.label = list(size=18))
ggsave(file = "../figures/phylomint_clsts.pdf", dpi=300, height=5, width=11)

plot.all = ggplot(phy.entero, aes(x=Classification, y= Distance, fill = Classification)) +
  geom_point_rast(colour = "darkgrey", size=0.5, position = position_jitterdodge(jitter.width = 0.5)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  scale_x_discrete(limits=c("Co-excluder", "Co-colonizer")) +
  scale_fill_manual(values = c("Co-excluder" = 'steelblue', "Co-colonizer" = "tomato")) +
  geom_signif(
    comparisons = list(c("Co-excluder", "Co-colonizer")),
    map_signif_level = FALSE) +
  ylab("Metabolic distance") +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12))
ggsave(plot.all, file = "../figures/phylomint_vs-entero.pdf", dpi=300, height=5, width=5)
