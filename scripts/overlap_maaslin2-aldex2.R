# load libraries
library(stringr)
library(tidyr)
library(data.table)
library(matrixStats)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
all.files = list.files("aldex2/", pattern=".tsv")
variables = str_match(all.files, "ALDEx2_\\s*(.*?)\\s*.tsv")

for (variable in variables[,2]) {
  cat(paste("Running analysis for", variable, "...\n"))
  aldex2.file = list.files(path = "aldex2/", pattern = paste(variable, ".tsv", sep=""), full.names=TRUE)
  maaslin2.file = list.files(path = "maaslin2/", pattern = paste(variable, ".tsv", sep=""), full.names=TRUE)
  genomes.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
  
  # filter data
  maaslin2.all = read.delim(maaslin2.file, stringsAsFactors = FALSE)
  meta.variable = str_split(string = variable, pattern = "_")[[1]][1]
  maaslin2.all = maaslin2.all[which(maaslin2.all$metadata == meta.variable),]
  maaslin2.all$qval = p.adjust(maaslin2.all$pval, method="fdr")
  maaslin2.sign = maaslin2.all[which(maaslin2.all$qval < 0.05),]
  aldex2.all = read.delim(aldex2.file, stringsAsFactors = FALSE)
  aldex2.sign = aldex2.all[which(aldex2.all$FDR < 0.05),]
  overlap = merge(maaslin2.sign, aldex2.sign, by="feature", all=FALSE)
  overlap = overlap[which(sign(overlap$coef) == sign(overlap[,colnames(aldex2.sign)[2]])),]
  
  # load taxonomy
  ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rownames(genomes.metadata) = genomes.metadata$Genome
  tax.df = unique(genomes.metadata[,c("Species_rep", "Lineage")])
  tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
  colnames(tax.df)[1] = "feature"
  overlap.tax = merge(overlap, tax.df, by="feature")
  species = str_split(string = meta.variable, pattern = "\\.")[[1]][1]
  if (species == "Entero") {
    overlap.tax = overlap.tax[which(overlap.tax$Family != "f__Enterobacteriaceae"),]
  } else if (species == "Ecoli") {
    overlap.tax = overlap.tax[which(overlap.tax$Genus != "g__Escherichia"),]
  } else if (species == "Kpneumo") {
    overlap.tax = overlap.tax[which(overlap.tax$Genus != "g__Klebsiella"),]
  }
  
  # report overlaps
  cat(paste("Found", nrow(overlap.tax), "overlapping species (Maaslin2 =", nrow(overlap.tax)/nrow(maaslin2.sign)*100, "; Aldex2 =", nrow(overlap.tax)/nrow(aldex2.sign)*100, ")\n"))
  
  # save table
  write.table(overlap.tax, file=paste0("overlap/", variable, ".tsv"), sep="\t", row.names=FALSE, quote=FALSE)
}