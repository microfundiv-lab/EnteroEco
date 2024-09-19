# load libraries
library(data.table)
library(stringr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggsignif)
library(CoDaSeq)
library(grid)
library(ggrastr)

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
metadata = read.delim("metadata/metagenomes_07-2024_samples.tsv", stringsAsFactors = FALSE, check.names = TRUE)
abund.data = fread("bwa/bwa_counts-filtered_batch-corr_samples.csv")
rownames(metadata) = metadata$Sample
abund.data.df = data.frame(abund.data[,-1], check.names = FALSE)
rownames(abund.data.df) = abund.data$Genome
abund.data.df = abund.data.df[,rownames(metadata)]

# load taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
genomes.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
rownames(genomes.metadata) = genomes.metadata$Genome
tax.df = unique(genomes.metadata[,c("Species_rep", "Lineage")])
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
rownames(tax.df) = tax.df$Species_rep
entero.species = as.vector(tax.df[which(tax.df$Family == "f__Enterobacteriaceae"),"Species_rep"])

# calculate Entero abundance
entero.row = colSums(abund.data.df[entero.species,])
abund.data.df2 = rbind(abund.data.df, entero.row)
rownames(abund.data.df2)[nrow(abund.data.df2)] = "Enterobacteriaceae"
abund.data.df2 = abund.data.df2[which(!rownames(abund.data.df2) %in% entero.species),]
abund.data.clr = codaSeq.clr(abund.data.df2 + 0.5, samples.by.row = FALSE)
metadata$Entero.abund = abund.data.clr["Enterobacteriaceae",]
metadata$Entero.presence = ifelse(metadata$Sample %in% colnames(abund.data.df2)[which(abund.data.df2["Enterobacteriaceae",] > 0)], "Yes", "No")

# format data
metadata = metadata[names(which(colSums(abund.data.df) > 0)),]
abund.data.df = abund.data.df[,rownames(metadata)]

# subsample
transpose.df = t(abund.data.df)
transpose.df = transpose.df[names(which(rowSums(transpose.df) >= 500000)),]
rarefy.df = rrarefy(transpose.df, 500000)

# calculate alpha diversity
metadata = metadata[intersect(rownames(rarefy.df), rownames(metadata)),]
metadata$Alpha.div = vegan::diversity(rarefy.df, index="shannon")
metadata$Richness = specnumber(rarefy.df)

# only health adults
metadata = metadata[which(metadata$Disease.name == "Healthy" & metadata$Age.group == "Adult"),]

# check correlation between diversity and entero abundance
corr = cor.test(metadata$Alpha.div, metadata$Entero.abund)
grob = grobTree(textGrob(paste("Pearson's R2 =", signif(corr$estimate**2, digits=3)), x=0.1,  y=0.8, hjust=0, gp=gpar(col="black", fontsize=12)))
entero.alpha = ggplot(metadata, aes(x=Entero.abund, y=Alpha.div)) +
  geom_point_rast(size=0.1) +
  geom_smooth(method="lm") +
  ylab("Alpha diversity (Shannon)") +
  xlab("Enterobacteriaceae abundance (CLR)") +
  theme_classic() +
  theme(axis.title.x = element_text(size=16)) +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  annotation_custom(grob)
ggsave("figures/alpha_diversity.pdf", dpi=300, height=4, width=9)
