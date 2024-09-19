# load libraries
library(reshape2)
library(vegan)
library(tidyr)
library(ggsignif)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# load data 
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/genofan/")
comb.abund = read.csv("../overlap/combined_coloniz-abund_filt.csv")
rownames(comb.abund) = comb.abund$Species_rep
comb.abund = comb.abund[which(comb.abund$Family != "Enterobacteriaceae"),]
comb.abund = comb.abund[which(comb.abund$N_genomes_c90 > 0),]
comb.abund = comb.abund[,-1]

# load kegg data
keggor.df = read.delim("kegg_orthologs.tsv", sep="", header = F)
keggor.df$V1 = gsub("_\\d+$", "", keggor.df$V1)
keggor.matrix = data.frame(acast(keggor.df, V2 ~ V1, length))
keggor.matrix = keggor.matrix[,rownames(comb.abund)]

# calculate shannon diversity
div.df = data.frame(diversity(t(keggor.matrix), index="shannon"))
div.df$V2 = comb.abund[rownames(div.df),"Classification"]
colnames(div.df) = c("Diversity", "Classification")

# load annotation data
annot.df = read.delim("annotation_coverage.tsv", header = F)
rownames(annot.df) = annot.df$V1
annot.df$coverage = annot.df$V4/annot.df$V3*100
annot.df = annot.df[rownames(comb.abund),]
annot.df$Classification = comb.abund[rownames(annot.df),"Classification"]

# plot boxplots
box.annot = ggplot(annot.df, aes(x=Classification, y=V4, fill=Classification)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_point(alpha=0.4, size=0.5, position = position_jitterdodge(jitter.width = 0.5)) +
  labs(y = "Number of annotated genes (KEGG)") +
  scale_fill_manual(values = c("Co-excluder" = 'steelblue', "Co-colonizer" = "tomato")) +
  scale_x_discrete(limits=c("Co-excluder", "Co-colonizer"), labels=c("Co-excluder", "Co-colonizer")) +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(
    comparisons = list(c("Co-excluder", "Co-colonizer")),
    map_signif_level = FALSE
  )

box.div = ggplot(div.df, aes(x=Classification, y=Diversity, fill=Classification)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_point(alpha=0.4, size=0.5, position = position_jitterdodge(jitter.width = 0.5)) +
  labs(y = "KO diversity (Shannon)") +
  scale_fill_manual(values = c("Co-excluder" = 'steelblue', "Co-colonizer" = "tomato")) +
  scale_x_discrete(limits=c("Co-excluder", "Co-colonizer"), labels=c("Co-excluder", "Co-colonizer")) +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(
    comparisons = list(c("Co-excluder", "Co-colonizer")),
    map_signif_level = FALSE
  )

# ggarrange
ggarrange(box.annot, box.div, common.legend=TRUE)
ggsave(file = "../figures/ko-wilcox_diversity-validation.pdf", dpi=300, height=4, width=7)