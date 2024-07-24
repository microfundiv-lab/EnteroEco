# load libraries
library(ape)
library(data.table)
library(vegan)
library(tidyr)
library(ggplot2)
library(ggsignif)

# input and output
cat("Loading data ...\n")
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
abund.data = fread("bwa/bwa_counts-filtered_batch-corr_samples.csv")
metadata = read.delim("metadata/metagenomes_07-2024_samples.tsv", check.names = TRUE)
rownames(metadata) = metadata$Sample
tax.df = unique(read.delim("metadata/species_uhgg_v1.2.tsv")[,c("Species_rep", "Lineage")])

# only healthy adults
metadata = metadata[which(metadata$Disease.name == "Healthy" & metadata$Age.group == "Adult"),]

# parse data
abund.data.df = data.frame(abund.data[,-1], check.names = FALSE)
rownames(abund.data.df) = abund.data$Genome
abund.data.df = abund.data.df[,metadata$Sample]
abund.data.df = abund.data.df[,which(colSums(abund.data.df) > 0)]
metadata = metadata[colnames(abund.data.df),]

# load taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
rownames(tax.df) = tax.df$Species_rep

# define entero presence/absence
entero.species = tax.df[which(tax.df$Family == "f__Enterobacteriaceae"),"Species_rep"]
entero.pos = names(which(colSums(abund.data.df[entero.species,]) > 0))
metadata$Entero.presence = ifelse(metadata$Sample %in% entero.pos, "Yes", "No")

# subsample
transpose.df = t(abund.data.df)
transpose.df = transpose.df[names(which(rowSums(transpose.df) >= 500000)),]
rarefy.df = rrarefy(transpose.df, 500000)

# generate distance matrix with CLR
cat("Calculating distance matrix ...\n")
dist.df.clr = coda.base::dist(rarefy.df+0.5, method="aitchison")
dist.df.clr = as.data.frame(as.matrix(dist.df.clr))

# calculate beta diversity per group
entero.pos.samples = intersect(rownames(metadata)[which(metadata$Entero.presence == "Yes")], rownames(dist.df.clr))
entero.neg.samples = intersect(rownames(metadata)[which(metadata$Entero.presence == "No")], rownames(dist.df.clr))
entero.pos.df = dist.df.clr[entero.pos.samples, entero.pos.samples]
entero.neg.df = dist.df.clr[entero.neg.samples, entero.neg.samples]

# prepare final df
entero.pos.values = as.data.frame(as.vector(as.matrix(entero.pos.df)))
colnames(entero.pos.values)[1] = "Distance"
entero.pos.values$Classification = "Presence"
entero.neg.values = as.data.frame(as.vector(as.matrix(entero.neg.df)))
colnames(entero.neg.values)[1] = "Distance"
entero.neg.values$Classification = "Absence"
beta.all = rbind(entero.pos.values, entero.neg.values)

# density plot
beta.dens = ggplot(beta.all, aes(x = Distance, fill = Classification)) +
  geom_density(alpha=0.5, colour="grey") +
  scale_fill_manual(values = c("Presence" = "tomato", "Absence" = 'steelblue'), name="Colonization status") +
  ylab("Density") +
  xlab("Pairwise Aitchison distance") +
  theme_classic() +
  theme(axis.title = element_text(size=14)) +
  theme(axis.text = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size=12))
ggsave(filename = "figures/beta-div_pairwise.pdf", width=8, height=4, dpi=300)