# load libraries
library(reshape2)
library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)

# load data 
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/gem")
comb.abund = read.csv("../overlap/combined_coloniz-abund_filt.csv")
rownames(comb.abund) = comb.abund$Species_rep
comb.abund = comb.abund[which(comb.abund$Family != "Enterobacteriaceae"),]
comb.abund = comb.abund[,-1]

# add genome type
genome.metadata = read.delim("../metadata/genomes_uhgg-v1.2.tsv")[,c("Genome", "Genome_type")]
rownames(genome.metadata) = genome.metadata$Genome

# load gem data
gem.df = read.delim("cobrapy-M3_results.tsv", header = F)
gem.df$V2 = make.names(gem.df$V2)

# function to plot comparison
compare_metab = function(df, type) {
  gem.df = df[which(df$V4 == type),]
  gem.matrix = data.frame(acast(gem.df, V1 ~ V2, length))
  gem.matrix = gem.matrix[rownames(comb.abund),]
  gem.matrix[gem.matrix > 1] = 1
  
  # calculate shannon diversity
  div.df = data.frame(specnumber(gem.matrix))
  div.df$V2 = comb.abund[rownames(div.df),"Classification"]
  colnames(div.df) = c("Diversity", "Classification")
  
  # plot boxplot
  box.plot = ggplot(div.df, aes(x=Classification, y=Diversity, fill=Classification)) +
    geom_boxplot(alpha=0.7, outlier.shape=NA) +
    geom_point(alpha=0.4, size=0.5, position = position_jitterdodge(jitter.width = 0.5)) +
    ylab(paste0("Number of metabolites (", tolower(type), ")")) +
    scale_fill_manual(values = c("Co-excluder" = 'steelblue', "Co-colonizer" = "tomato")) +
    guides(fill="none") +
    scale_x_discrete(limits=c("Co-excluder", "Co-colonizer")) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size=14)) +
    theme(axis.text.x = element_text(size=12)) +
    theme(axis.text.y = element_text(size=12)) +
    geom_signif(
      comparisons = list(c("Co-excluder", "Co-colonizer")),
      map_signif_level = FALSE
    )
  return(box.plot)
}

# plot
box.up = compare_metab(gem.df, "Uptake")
box.secr = compare_metab(gem.df, "Secretion")
ggarrange(box.up, box.secr, ncol=1)
ggsave(file = "../figures/gem_richness.pdf", dpi=300, height=10, width=4)
