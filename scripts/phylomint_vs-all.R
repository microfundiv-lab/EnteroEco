# load libraries
library(data.table)
library(tidyverse)
library(ggsignif)
library(ggrastr)
library(ggpubr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/gem")
phylomint.res = read.table("phylomint_gapless.tsv", fill=TRUE, header=TRUE)
phylomint.res$identifier = apply(phylomint.res[, c("A", "B")], 1, function(x) paste(sort(x), collapse="_"))
phylomint.res = phylomint.res[!duplicated(phylomint.res$identifier), ]
phylomint.res$identifier = NULL
metadata = read.csv("../overlap/combined_coloniz-abund_filt.csv")
rownames(metadata) = metadata$Species_rep
metadata = metadata[which(metadata$Family != "Enterobacteriaceae"),]

# get list of species
co.excluders = rownames(metadata)[which(metadata$Classification == "Co-excluder")]
co.colonizers = rownames(metadata)[which(metadata$Classification == "Co-colonizer")]
entero.species = read.delim("../metadata/species_entero.txt", header = F, stringsAsFactors = F)

# add taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
species.data = read.delim("../metadata/species_uhgg_v1.2.tsv")
species.data = unique(species.data[,c("Species_rep", "Lineage")])
species.data = separate(data = species.data, col = Lineage, sep = ";", into = ranks)
rownames(species.data) = species.data$Species_rep
species.data = data.frame(lapply(species.data, function(x) { gsub(".__", "", x)}))
rownames(species.data) = species.data$Species_rep

# function for comparison
phylo_compare = function(target) {
  phy.df = phylomint.res[which((phylomint.res$A %in% metadata$Species_rep & phylomint.res$B %in% target) | (phylomint.res$B %in% metadata$Species_rep & phylomint.res$A %in% target)),]
  phy.df$Species_rep = ifelse(phy.df$A %in% target, phy.df$B, phy.df$A)
  phy.df$Target = ifelse(phy.df$A == phy.df$Species_rep, phy.df$B, phy.df$A)
  phy.df = phy.df[which(phy.df$Species_rep != phy.df$Target),]
  return(phy.df)
}
phy.colon = phylo_compare(co.colonizers)
phy.colon$Comparison = "vs. Co-colonizers"
phy.coexcl = phylo_compare(co.excluders)
phy.coexcl$Comparison = "vs. Co-excluders"
phy.fi = rbind(phy.colon, phy.coexcl)

# calculate distances
phy.fi$Distance = 1-(as.numeric(phy.fi$Competition) - as.numeric(phy.fi$Complementarity))
phy.fi = merge(phy.fi, metadata[c("Species_rep", "Classification", "Order", "Family", "Genus")], by = "Species_rep")

# plot and save
plot.split = ggplot(phy.fi, aes(x=Classification, y= Distance, fill = Classification)) +
  geom_point_rast(alpha=0.1, size=0.5, position = position_jitterdodge(jitter.width = 0.5)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~ Comparison, ncol=3) +
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
ggsave(file = "../figures/phylomint_vs-all.pdf", dpi=300, height=5, width=7)
