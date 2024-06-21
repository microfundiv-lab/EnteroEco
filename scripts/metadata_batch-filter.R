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
metadata = read.delim("metadata/metagenomes_11-2023_samples.tsv", check.names = TRUE)
rownames(metadata) = metadata$Sample
tax.df = unique(read.delim("metadata/species_uhgg_v1.2.tsv")[,c("Species_rep", "Lineage")])

# filter metadata
metadata = metadata[which(!is.na(metadata$Disease.name) & !is.na(metadata$Age.group) & !is.na(metadata$Country)),]

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
studies = rownames(batch.prop)
studies = studies[which(studies != "SRP001634")]

# clean metadata
metadata = metadata[which(metadata$Study %in% studies & !metadata$Sample %in% c("SRR3340629", "SRR3340631")),]
write.table(metadata, file="metadata/metagenomes_11-2023.tsv", sep="\t", row.names=FALSE, quote=FALSE)