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
process_phylomint = function(input) {
  diet = gsub(".tsv", "", gsub(".*_", "", input))
  phylomint.res = read.table(input, fill=TRUE, header=TRUE)
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
  phy.entero$Medium = diet
  return(phy.entero)
}

# run function for all diets
diets = list.files(".", "phylomint*")
phylomint.list = lapply(diets, function(x) {
  cat("Generating table for", x, "...\n")
  phy.df = process_phylomint(x)
  return(phy.df)
})
phy.combined = as.data.frame(rbindlist(phylomint.list))

# format diets
phy.combined = phy.combined[which(phy.combined$Medium != "gapless"),]
phy.combined$Medium = gsub("average", "EU average", phy.combined$Medium)
phy.combined$Medium = gsub("dach", "DACH", phy.combined$Medium)
phy.combined$Medium = gsub("fiber", "High fiber", phy.combined$Medium)
phy.combined$Medium = gsub("free", "Gluten free", phy.combined$Medium)
phy.combined$Medium = gsub("lowcarb", "High fat, low carb", phy.combined$Medium)
phy.combined$Medium = gsub("medite", "Mediterranean", phy.combined$Medium)
phy.combined$Medium = gsub("protein", "High protein", phy.combined$Medium)
phy.combined$Medium = gsub("t2d", "T2D", phy.combined$Medium)
phy.combined$Medium = gsub("unhealthy", "Unhealthy", phy.combined$Medium)
phy.combined$Medium = gsub("vegan", "Vegan", phy.combined$Medium)
phy.combined$Medium = gsub("veget", "Vegetarian", phy.combined$Medium)

# plot results
plot.all = ggplot(phy.combined, aes(x=Classification, y= Distance, fill = Classification)) +
  geom_point_rast(colour = "darkgrey", size=0.5, position = position_jitterdodge(jitter.width = 0.5)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~ Medium) +
  theme_classic() +
  scale_x_discrete(limits=c("Co-excluder", "Co-colonizer")) +
  scale_fill_manual(values = c("Co-excluder" = 'steelblue', "Co-colonizer" = "tomato")) +
  ylab("Metabolic distance") +
  theme(strip.text = element_text(size=16)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=14))
ggsave(filename="../figures/phylomint_diets.pdf", dpi=300, height=12, width=18)
