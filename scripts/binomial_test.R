# load libraries
library(data.table)
library(tidyr)
library(ggrastr)
library(ggplot2)

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
abund.data = fread("bwa/bwa_counts-filtered_batch-corr_samples.csv")
abund.data.df = data.frame(abund.data[,-1], check.names = FALSE)
rownames(abund.data.df) = abund.data$Genome

# filter metadata
metadata = read.delim("metadata/metagenomes_11-2023_samples.tsv", stringsAsFactors = FALSE, check.names = TRUE)
rownames(metadata) = metadata$Sample
abund.data.df = abund.data.df[,metadata$Sample]

# load taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
genomes.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
rownames(genomes.metadata) = genomes.metadata$Genome
tax.df = unique(genomes.metadata[,c("Species_rep", "Lineage")])
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
rownames(tax.df) = tax.df$Species_rep

# define species presence/absence
ecoli.species = "GUT_GENOME144544" # E.coli
ecoli.pos = names(which(colSums(abund.data.df[ecoli.species,]) > 0))
kpneumo.species = "GUT_GENOME147598"
kpneumo.pos = names(which(colSums(abund.data.df[kpneumo.species,]) > 0))
enrogg.species = "GUT_GENOME143746"
enrogg.pos = names(which(colSums(abund.data.df[enrogg.species,]) > 0))
metadata$Ecoli.presence = ifelse(metadata$Sample %in% ecoli.pos, "Yes", "No")
metadata$Kpneumo.presence = ifelse(metadata$Sample %in% kpneumo.pos, "Yes", "No")
metadata$Enrogg.presence = ifelse(metadata$Sample %in% enrogg.pos, "Yes", "No")

# double co-cocolonization
ecoli_freq = length(which(metadata$Ecoli.presence == "Yes"))/nrow(metadata)
kpneumo_freq = length(which(metadata$Kpneumo.presence == "Yes"))/nrow(metadata)
expprob = ecoli_freq*kpneumo_freq
obsfreq = length(which(metadata$Ecoli.presence == "Yes" & metadata$Kpneumo.presence == "Yes"))
stat = binom.test(obsfreq, nrow(metadata), expprob)

# triple co-colonization
enrogg_freq = length(which(metadata$Enrogg.presence == "Yes"))/nrow(metadata)
expprob_3 = ecoli_freq*kpneumo_freq*enrogg_freq
obsfreq_3 = length(which(metadata$Ecoli.presence == "Yes" & metadata$Kpneumo.presence == "Yes" & metadata$Enrogg.presence == "Yes"))
stat_3 = binom.test(obsfreq_3, nrow(metadata), expprob_3)
