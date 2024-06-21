# load libraries
library(reshape2)
library(vegan)
library(tidyr)
library(ggsignif)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# load data 
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/genofan/")
comb.abund = read.csv("../overlap/combined_coloniz-abund_filt.csv")
rownames(comb.abund) = comb.abund$Species_rep
comb.abund = comb.abund[which(comb.abund$Family != "Enterobacteriaceae"),]
comb.abund = comb.abund[,-1]

# parse kegg results
ko.classes.ori = read.delim("KO_Orthology_ko00001.txt", header=FALSE)
ko.classes = separate(ko.classes.ori, V4, "V5", sep=" ")
ko.desc = cbind(ko.classes.ori, ko.classes)
ko.desc = unique(ko.desc[,c("V5", "V4")])
rownames(ko.desc) = ko.desc$V5

# load kegg data
keggor.df = read.delim("kegg_orthologs.tsv", sep="", header = F)
keggor.df$V1 = gsub("_\\d+$", "", keggor.df$V1)
keggor.matrix = data.frame(acast(keggor.df, V2 ~ V1, length))
keggor.matrix = keggor.matrix[,rownames(comb.abund)]

# calculate shannon diversity
div.df = data.frame(diversity(t(keggor.matrix)))
div.df$V2 = comb.abund[rownames(div.df),"Classification"]
colnames(div.df) = c("Diversity", "Classification")

# plot boxplot
box.plot = ggplot(div.df, aes(x=Classification, y=Diversity, fill=Classification)) +
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
ggsave(file = "../figures/ko-wilcox_diversity.pdf", dpi=300, height=5, width=5)