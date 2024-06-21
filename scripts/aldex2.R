#!/usr/bin/env Rscript

# load libraries
library(data.table)
library(stringr)
library(tidyr)
library(ALDEx2)

# arguments
args = commandArgs(trailingOnly=TRUE)

# input and output
metadata = read.delim("/rds/project/rds-aFEMMKDjWlo/aalmeida/project_ecoevo/metadata/metagenomes_11-2023_samples.tsv", stringsAsFactors = FALSE, check.names = TRUE)
abund.data = fread("/rds/project/rds-aFEMMKDjWlo/aalmeida/project_ecoevo/bwa/bwa_counts-filtered_samples.csv")
genomes.metadata = read.delim("/rds/project/rds-aFEMMKDjWlo/rfs_data/metadata/species_uhgg_v1.2.tsv", stringsAsFactors = FALSE)
variable = args[1] # Entero.presence, Ecoli.presence, Kpneumo.presence
meta.filter = args[2] # All, Healthy-Adult
outfile = paste("ALDEx2_", variable, "_", meta.filter, ".tsv", sep="")

# format data
rownames(metadata) = metadata$Sample
abund.data.df = data.frame(abund.data[,-1], check.names = FALSE)
rownames(abund.data.df) = abund.data$Genome
abund.data.df = abund.data.df[,metadata$Sample]
abund.data.df = abund.data.df[,which(colSums(abund.data.df) > 0)]
metadata = metadata[colnames(abund.data.df),]

# load taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(genomes.metadata) = genomes.metadata$Genome
tax.df = separate(data = genomes.metadata, col = Lineage, sep = ";", into = ranks)

# define entero presence/absence
entero.species = tax.df[which(tax.df$Family == "f__Enterobacteriaceae"),"Species_rep"]
entero.pos = names(which(colSums(abund.data.df[entero.species,]) > 0))
metadata$Entero.presence = ifelse(metadata$Sample %in% entero.pos, "Yes", "No")

# calculate relative abundance
abund.data.df = abund.data.df[tax.df$Species_rep,]
bwa.relab = abund.data.df/(tax.df$Length/1000)
for(col in colnames(bwa.relab)) {
  bwa.relab[,col] = bwa.relab[,col] / sum(bwa.relab[,col]) * 100
}
metadata$Ecoli.abund = colSums(bwa.relab["GUT_GENOME144544",])
metadata$Ecoli.presence = ifelse(metadata$Ecoli.abund > 0, "Yes", "No")
metadata$Kpneumo.abund = colSums(bwa.relab["GUT_GENOME147598",])
metadata$Kpneumo.presence = ifelse(metadata$Kpneumo.abund > 0, "Yes", "No")
res_column = paste(variable, "Yes", sep="")

# filter metadata based on arguments
metadata$Disease.name = ifelse(metadata$Health.state == "Healthy", "Healthy", metadata$Disease.name)
metadata$Age.group = relevel(factor(metadata$Age.group), "Adult")
metadata$Disease.name = relevel(factor(metadata$Disease.name), "Healthy")
metadata$Continent = relevel(factor(metadata$Continent), "Europe")

if (meta.filter == "All") {
  metadata = metadata[which(!is.na(metadata$Disease.name) & !is.na(metadata$Age.group) & !is.na(metadata$Continent)),]
  pos.studies = unique(metadata[which(metadata[,variable] == "Yes"),"Study"])
  neg.studies = unique(metadata[which(metadata[,variable] == "No"),"Study"])
  studies = intersect(pos.studies, neg.studies)
  metadata = metadata[which(metadata$Study %in% studies),]
  mm = model.matrix(formula(paste("~", variable, "+ Disease.name + Age.group + Continent + Read.count + Study")), metadata)
} else if (meta.filter == "Healthy-Adult") {
  metadata = metadata[which(metadata$Age.group == "Adult" & metadata$Disease.name == "Healthy" & !is.na(metadata$Continent)),]
  pos.studies = unique(metadata[which(metadata[,variable] == "Yes"),"Study"])
  neg.studies = unique(metadata[which(metadata[,variable] == "No"),"Study"])
  studies = intersect(pos.studies, neg.studies)
  metadata = metadata[which(metadata$Study %in% studies),]
  mm = model.matrix(formula(paste("~", variable, "+ Continent + Read.count + Study")), metadata)
}

abund.data.df = abund.data.df[,rownames(mm)]

# prevalence filter
thresh = 0.01 # 1% prevalence
filtered_genomes = names(which(rowSums(abund.data.df > 0)/ncol(abund.data.df) > thresh))
abund.data.filtered = abund.data.df[filtered_genomes,]

# run ALDEx2
cat("Running ALDEx2 with", ncol(abund.data.filtered), "samples and", nrow(abund.data.filtered), "species ...")
clr.values = aldex.clr(abund.data.filtered, mm, mc.samples = 32, denom = "all", useMC = TRUE) # genomes in rows, samples in columns
glm.test = aldex.glm(clr.values, verbose = TRUE)
all.results = glm.test[,c(paste(res_column,":Est", sep=""), paste(res_column, ":pval", sep=""))]
all.results$FDR = p.adjust(all.results[,2], method="fdr")

# save table
all.results.fi = cbind(rownames(all.results), all.results)
colnames(all.results.fi)[1] = "feature"
write.table(all.results.fi, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

