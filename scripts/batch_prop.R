# load libraries
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(reshape2)
library(data.table)
library(grid)

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
abund.data = fread("bwa/bwa_counts-filtered_samples.csv")
metadata = read.delim("metadata/metagenomes_07-2024_samples.tsv", check.names = TRUE)
rownames(metadata) = metadata$Sample
tax.df = unique(read.delim("metadata/species_uhgg_v1.2.tsv")[,c("Species_rep", "Lineage")])

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

# split dataset into 
batch.prop = as.data.frame.matrix(table(metadata[,c("Study", "Entero.presence")]))
batch.prop$Diff = batch.prop$Yes-batch.prop$No
batch.prop = batch.prop[which(batch.prop$Yes > 0 & batch.prop$No > 0),]

# analyse data
correlation = cor.test(log10(batch.prop$No+1), log10(batch.prop$Yes+1), alternative="two.sided")
grob = grobTree(textGrob(paste("Pearson's R2 =", signif(correlation$estimate**2, digits=3), "\nP =", 
                                          signif(correlation$p.value, digits=3)), x=0.2,  y=0.8, hjust=0, gp=gpar(col="black", fontsize=14)))
scatter.plot = ggplot(batch.prop, aes(x=log10(Yes+1), y=log10(No+1))) +
  geom_point() +
  geom_smooth(method="lm") +
  annotation_custom(grob) +
  ylab(bquote("Samples w/ Enterobacteriaceae (log"[10]*")")) +
  xlab(bquote("Samples w/o Enterobacteriaceae (log"[10]*")")) +
  theme_classic() +
  theme(axis.title.x = element_text(size=16)) +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14))
ggsave("figures/batch_props.pdf", height=4, width=12)
