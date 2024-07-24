# load libraries
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(scales)
library(stringr)
library(data.table)
library(reshape2)

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data")
overlap_dir = "overlap/"

# load and process taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
species.metadata = read.delim("metadata/species_uhgg_v1.2.tsv", stringsAsFactors = FALSE)
rownames(species.metadata) = species.metadata$Genome
tax.df = unique(species.metadata[,c("Species_rep", "Lineage")])
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
tax.df$Taxon = NA
for (n in 1:nrow(tax.df)){
  for (r in ranks){
    if(!grepl("__$", tax.df[n,r])){
      tax.df[n,"Taxon"] = tax.df[n,r] }}}

tax.df = data.frame(lapply(tax.df, function(x) { gsub(".__", "", x)}))
rownames(tax.df) = tax.df$Species_rep

# load overlap results
input.files =  list.files(path=overlap_dir, pattern = ".tsv", full.names = TRUE)
analyses.list = lapply(input.files, function(i) {
  analysis = strsplit(basename(i), split="\\.tsv")[[1]][1]
  target = strsplit(analysis, split="\\.")[[1]][1]
  df = read.delim(i); df$analysis = analysis
  colnames(df)[which(colnames(df) == paste0(target,".presenceYes.Est"))] = "ALDEx2:Est"
  colnames(df)[which(colnames(df) == paste0(target,".presenceYes.pval"))] = "ALDEx2:pval"
  return(df)
})

# combine all results into a dataframe
all.results = as.data.frame(rbindlist(analyses.list))

# rescale data
all.results$coef = all.results$coef/max(abs(all.results$coef))
all.results[,"ALDEx2:Est"] = all.results[,"ALDEx2:Est"]/max(abs(all.results[,"ALDEx2:Est"]))
all.results$overlap_mean = rowMeans(all.results[,c("coef", "ALDEx2:Est")])

# reshape data
overlap.df = as.data.frame(acast(feature ~ analysis, value.var="overlap_mean", data=all.results))
overlap.df$Genomes = rownames(overlap.df)

# function to read all fastspar result files
process_fastspar = function(corr_df, pval_df, metavariable, taxon=NULL) {
  d.corr_full = read.table(corr_df, sep='\t', header=TRUE, comment.char='', row.names=1, check.names=FALSE)
  d.pval_full = read.table(pval_df, sep='\t', header=TRUE, comment.char='', row.names=1, check.names=FALSE)
  
  # Mask upper triangle with NAs then melt, excluding upper triangle values (no repeating value)
  d.corr_full[upper.tri(d.corr_full, diag=TRUE)] <- NA
  d.pval_full[upper.tri(d.pval_full, diag=TRUE)] <- NA
  d.corr = reshape2::melt(as.matrix(d.corr_full), na.rm=TRUE, varnames=c('otu_1', 'otu_2'), value.name='correlation')
  d.pval = reshape2::melt(as.matrix(d.pval_full), na.rm=TRUE, varnames=c('otu_1', 'otu_2'), value.name='pvalue')
  
  # Merge correlations and pvalues
  d = merge(d.corr, d.pval)
  d$FDR = p.adjust(d$pvalue, method="BH")
  d = d[which(d$FDR < 0.05),]
  d$otu_1 = as.character(d$otu_1)
  d$otu_2 = as.character(d$otu_2)
  
  if (is.null(taxon)) {
    # Split into Ecoli and Kpneumo
    ecoli.df = d[which(d$otu_1 == "GUT_GENOME144544" | d$otu_2 == "GUT_GENOME144544"),]
    ecoli.df$neighbour = ifelse(ecoli.df$otu_1 == "GUT_GENOME144544", ecoli.df$otu_2, ecoli.df$otu_1)
    kpneumo.df = d[which(d$otu_1 == "GUT_GENOME147598" | d$otu_2 == "GUT_GENOME147598"),]
    kpneumo.df$neighbour = ifelse(kpneumo.df$otu_1 == "GUT_GENOME144544", kpneumo.df$otu_2, kpneumo.df$otu_1)
    
    # combine dataframes
    sign.species = unique(c(ecoli.df$neighbour, kpneumo.df$neighbour))
    combined.df = data.frame(matrix(nrow=length(sign.species), ncol=3))
    rownames(combined.df) = sign.species
    combined.df[,1] = sign.species
    colnames(combined.df) = c("Genomes", "Ecoli.fastspar", "Kpneumo.fastspar")
    
    # get correlation values
    combined.df$Ecoli.fastspar = ecoli.df[match(rownames(combined.df), ecoli.df$neighbour),"correlation"]
    combined.df$Kpneumo.fastspar = kpneumo.df[match(rownames(combined.df), kpneumo.df$neighbour),"correlation"]
    colnames(combined.df) = c("Genomes", paste("Ecoli.fastspar", metavariable, sep="_"), paste("Kpneumo.fastspar", metavariable, sep="_"))
    fastspar.df = combined.df
  } else {
    # select Enterobacteriaceae
    entero.df = d[which(d$otu_1 == "Enterobacteriaceae" | d$otu_2 == "Enterobacteriaceae"),]
    entero.df$neighbour = ifelse(entero.df$otu_1 == "Enterobacteriaceae", entero.df$otu_2, entero.df$otu_1)
    
    # combine into final dataframe
    fastspar.df = data.frame(matrix(ncol=1, nrow=length(unique(entero.df$neighbour))))
    rownames(fastspar.df) = unique(entero.df$neighbour)
    colnames(fastspar.df) = paste("Entero.fastspar", metavariable, sep="_")
    fastspar.df[match(entero.df$neighbour, rownames(fastspar.df)),1] = entero.df$correlation
    fastspar.df$Genomes = rownames(fastspar.df)
  }
  return(fastspar.df)
}

fastspar.healthy = process_fastspar("fastspar/healthy-adult_species_corr.tsv", "fastspar/healthy-adult_species_pvalues.tsv", "Healthy-Adult")
fastspar.highqual = process_fastspar("fastspar/all_species_corr.tsv", "fastspar/all_species_pvalues.tsv", "All")
fastspar.entero.healthy = process_fastspar("fastspar/healthy-adult_entero_corr.tsv", "fastspar/healthy-adult_entero_pvalues.tsv", "Healthy-Adult", "Entero")
fastspar.entero.highqual = process_fastspar("fastspar/all_entero_corr.tsv", "fastspar/all_entero_pvalues.tsv", "All", "Entero")

# Mmrge all fastspar dataframes
df_list = list(fastspar.highqual, fastspar.healthy, fastspar.entero.highqual, fastspar.entero.healthy)
abund.df = Reduce(function(x, y) merge(x, y, by="Genomes", all=TRUE), df_list)
abund.df[is.na(abund.df)] = 0
abund.df[,-1] = abund.df[,-1]/max(abs(abund.df[,-1]))

# merge colonization and abundance
combined.df = merge(overlap.df, abund.df, by="Genomes")
rownames(combined.df) = combined.df$Genomes
combined.df[is.na(combined.df)] = 0
colnames(combined.df)[1] = "Species_rep"

# filter for coloniz and abund
taxa = c("Ecoli", "Kpneumo", "Entero")
meta = c("All", "Healthy-Adult")
selected = c()
for (t in taxa){
  for (m in meta){
   pos.genomes = names(which(rowSums(combined.df[,grep(paste0(t,".*",m), colnames(combined.df))] > 0) == 2))
   neg.genomes = names(which(rowSums(combined.df[,grep(paste0(t,".*",m), colnames(combined.df))] < 0) == 2))  
   selected = c(selected, pos.genomes, neg.genomes)  
  }
}
selected = unique(selected)
final.df = combined.df[selected,]

# reformat columns
order.cols = c("Entero.presence_All", "Entero.fastspar_All", "Ecoli.presence_All", "Ecoli.fastspar_All", "Kpneumo.presence_All", "Kpneumo.fastspar_All",
               "Entero.presence_Healthy-Adult", "Entero.fastspar_Healthy-Adult", "Ecoli.presence_Healthy-Adult", "Ecoli.fastspar_Healthy-Adult", "Kpneumo.presence_Healthy-Adult", "Kpneumo.fastspar_Healthy-Adult")

final.df = final.df[,c("Species_rep", order.cols)]
colnames(final.df) = gsub(".presence", ".colonization", colnames(final.df))
colnames(final.df) = gsub(".fastspar", ".abundance", colnames(final.df))

# calculate summary scores
final.df$`Co-excluder`= rowSums(final.df[,2:13] < 0)
final.df$`Co-colonizer` = rowSums(final.df[,2:13] > 0)
final.df$ES = rowSums(final.df[,2:13])/(rowSums(final.df[,2:13] == 0)+1)

# save final file
genome.metadata = read.delim("metadata/genomes_uhgg-v1.2.tsv")
rownames(genome.metadata) = genome.metadata$Genome
genome.c90 = genome.metadata[which(genome.metadata$Completeness > 90),]
genome.c90 = data.frame(table(genome.c90$Species_rep))
colnames(genome.c90) = c("Species_rep", "N_genomes_c90")
metadata = merge(final.df, tax.df, by="Species_rep")
metadata$Classification = ifelse(metadata$ES > 0, "Co-colonizer", "Co-excluder")
rownames(metadata) = metadata$Row.names
metadata.c90 = merge(metadata, genome.c90, by="Species_rep", all.x = TRUE)
metadata.c90$N_genomes_c90 = ifelse(is.na(metadata.c90$N_genomes_c90), 0, metadata.c90$N_genomes_c90)
metadata.c90 = metadata.c90[,which(colnames(metadata.c90) != "Row.names")]

# save files
exclude = rownames(metadata.c90)[which(metadata.c90$`Co-excluder` > 0 & metadata.c90$`Co-colonizer` > 0)]
metadata.c90.filt = metadata.c90[!rownames(metadata.c90) %in% exclude,]
write.table(metadata.c90, file = paste0(overlap_dir, "/combined_coloniz-abund_all.csv"), row.names=FALSE, quote=FALSE, sep=",")
write.table(metadata.c90.filt, file = paste0(overlap_dir, "/combined_coloniz-abund_filt.csv"), row.names=FALSE, quote=FALSE, sep=",")

# save candidate files (species)
metadata.c90.only = metadata.c90.filt[which(metadata.c90.filt$N_genomes_c90 >= 10),]
write.table(genome.metadata[metadata.c90.filt[,1],"FTP_download"], file=paste0(overlap_dir, "/candidate_species_paths.txt"), row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(metadata.c90.filt[,1], file=paste0(overlap_dir, "/candidate_species.txt"), row.names=FALSE, quote=FALSE, col.names=FALSE)

# save candidate files (strains)
entero.health.sig = which(metadata.c90.only[,"Entero.colonization_Healthy-Adult"] != 0 & metadata.c90.only[,"Entero.abundance_Healthy-Adult"] != 0)
metadata.c90.strains = metadata.c90.only[entero.health.sig,]
write.table(metadata.c90.strains[,1], file=paste0(overlap_dir, "/selected_strains.txt"), row.names=FALSE, quote=FALSE, col.names=FALSE)