# load libraries
library(reshape2)
library(ggplot2)
library(matrixStats)
library(tidyr)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggthemes)
library(ggtext)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
metadata = read.delim("metadata/metagenomes_07-2024_samples.tsv", stringsAsFactors = FALSE, check.names = TRUE)
metadata = metadata[which(!is.na(metadata$Disease.name) & !is.na(metadata$Age.group) & !is.na(metadata$Continent)),]
rownames(metadata) = metadata$Sample

# rename diseases
metadata$Disease.name = gsub("CRC", "Colorectal cancer", metadata$Disease.name)
metadata$Disease.name = gsub("CD", "Crohn's disease", metadata$Disease.name)
metadata$Disease.name = gsub("UC", "Ulcerative colitis", metadata$Disease.name)
metadata$Disease.name = gsub("T1D", "Type 1 diabetes", metadata$Disease.name)
metadata$Disease.name = gsub("T2D", "Type 2 diabetes", metadata$Disease.name)
metadata$Disease.name = gsub("T2D", "Type 2 diabetes", metadata$Disease.name)
metadata$Disease.name = gsub("AS", "Ankylosing spondylitis", metadata$Disease.name)
metadata$Disease.name = gsub("RA", "Rheumatoid arthritis", metadata$Disease.name)
metadata$Disease.name = gsub("Parkinson", "Parkinson's disease", metadata$Disease.name)
metadata$Disease.name = gsub("GVHD", "Graft-versus-host disease", metadata$Disease.name)
metadata$Disease.name = gsub("MS", "Multiple sclerosis", metadata$Disease.name)
metadata$Disease.name = gsub("ME", "Myalgic encephalomyelitis", metadata$Disease.name)
metadata$Disease.name = gsub("BoneDisease", "Bone disease", metadata$Disease.name)

# load genome metadata
genomes.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
rownames(genomes.metadata) = genomes.metadata$Genome
genome.lengths = genomes.metadata[unique(as.vector(genomes.metadata$Species_rep)),c("Genome", "Length")]
gtdb.tax = unique(genomes.metadata[,c("Species_rep", "Lineage")])
rownames(gtdb.tax) = gtdb.tax$Species_rep
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
gtdb.tax = separate(data = gtdb.tax, col = Lineage, sep = ";", into = ranks)
species.data = gtdb.tax[,-1]
family.spcies.data = species.data[,c("Family","Species")]
family.spcies.data$Genome = rownames(family.spcies.data)

# load count data
abund.data = as.data.frame(fread("bwa/bwa_counts-filtered_samples.csv"))
rownames(abund.data) = abund.data$Genome
abund.data = abund.data[which(abund.data$Genome %in% genome.lengths$Genome),]
abund.data.df = abund.data[,rownames(metadata)]
metadata = metadata[colnames(abund.data.df),]
genome.lengths = genome.lengths[rownames(abund.data.df),]

# select species
entero.species = family.spcies.data[which(family.spcies.data$Family == "f__Enterobacteriaceae"),"Genome"]
overall.prev = data.frame(sort(rowSums(abund.data.df > 0), decreasing=TRUE))
overall.prev = overall.prev[which(overall.prev[,1] > 0), ,drop=FALSE]
detected.entero = rownames(overall.prev)[rownames(overall.prev) %in% entero.species]

# function to estimate prevalence
prevalence_by_metadata = function(prevalence_data, metadata.df, variable) {
  df.variable = data.frame(matrix(ncol=length(unique(metadata.df[,variable])), nrow=nrow(prevalence_data)))
  rownames(df.variable) = rownames(prevalence_data)
  colnames(df.variable) = unique(metadata.df[,variable])
  
  for (value in colnames(df.variable)){
    selected.samples = metadata.df[which(metadata.df[,variable] == value),"Sample"]
    df.filtered = prevalence_data[,which(colnames(prevalence_data) %in% selected.samples)]
    df.prevalence = data.frame(rowSums(df.filtered > 0)/ncol(df.filtered)*100)
    df.variable[rownames(df.prevalence),value] = df.prevalence
    df.variable$Genome = rownames(df.variable)
  }
  return(df.variable)
}

# calculate prevalence per continent, age and disease
df.continent = prevalence_by_metadata(prevalence_data = abund.data.df, metadata.df = metadata, variable = "Continent")
df.age = prevalence_by_metadata(prevalence_data = abund.data.df, metadata.df = metadata, variable = "Age.group")
df.health = prevalence_by_metadata(prevalence_data = abund.data.df, metadata.df = metadata, variable = "Disease.name")

# function to estimate relative abundance
relab_by_metadata <- function(relab_data, metadata, genome_lengths, variable) {
  relab_variable <- data.frame(matrix(ncol = length(unique(metadata[, variable])), nrow = nrow(relab_data)))
  rownames(relab_variable) <- rownames(relab_data)
  colnames(relab_variable) <- unique(metadata[, variable])
  
  for (value in colnames(relab_variable)) {
    selected_samples <- metadata[which(metadata[, variable] == value), "Sample"]
    relab_filtered <- relab_data[, selected_samples]
    relab_filtered = relab_filtered[rownames(genome_lengths),]
    bwa.relab = relab_filtered/(genome_lengths$Length/1000)
    
    for(col in colnames(bwa.relab)) {
      bwa.relab[,col] = bwa.relab[,col] / sum(bwa.relab[,col]) * 100
    }
    # Calculate the average abundance of each genome
    bwa.relab[bwa.relab == 0] = NA
    median <- rowMedians(as.matrix(bwa.relab), na.rm=TRUE)
    relab_variable[,value] = median
  }
  relab_variable$Genome = rownames(relab_variable)
  return(relab_variable)
}

# prepara dataframe, calculate rel_abundance of genome in each condition and category
relab.age = relab_by_metadata(relab_data = abund.data.df, metadata = metadata, genome_lengths = genome.lengths, variable = "Age.group")
relab.health = relab_by_metadata(relab_data = abund.data.df, metadata = metadata, genome_lengths = genome.lengths, variable = "Disease.name")
relab.continent = relab_by_metadata(relab_data = abund.data.df, metadata = metadata, genome_lengths = genome.lengths, variable = "Continent")

# select Enterobacteriaceae
generate_final_dataframe = function(condition.value, df.input, relab.input){
  # Melt the matrices into long format
  df.condition.melted = reshape2::melt(df.input, variable.name = "condition.value", value.name = "Prevalence")
  relab.condition.melted = reshape2::melt(relab.input, variable.name = "condition.value", value.name = "Rel_abund")
  
  # Merge the melted data frames
  combined_data <- merge(df.condition.melted, relab.condition.melted, by = c("Genome", "condition.value"))
  colnames(combined_data)[1] = "Genome"
  combined_data = merge(combined_data, family.spcies.data[,c("Species","Genome")], by = "Genome")
  
  # remove s__ prefix
  combined_data$Species = gsub("s__","",combined_data$Species)
  
  # final output
  combined_data$Metadata = condition.value
  return(combined_data)
}

final_age = generate_final_dataframe("Age.group", df.age, relab.age)
final_health = generate_final_dataframe("Disease.name", df.health, relab.health)
final_continent =generate_final_dataframe("Continent", df.continent, relab.continent)
final_df = rbind(final_age, final_health, final_continent)

# rename metadata
final_df$Metadata = gsub("Age.group", "Age group", final_df$Metadata)
final_df$Metadata = gsub("Disease.name", "Health state", final_df$Metadata)

# order alphabetically
final_df$condition.value = factor(final_df$condition.value, levels=sort(as.vector(unique(final_df$condition.value))))

# create balloon plot
final_df = final_df[which(final_df$Genome %in% c("GUT_GENOME144544", "GUT_GENOME147598","GUT_GENOME143760", "GUT_GENOME000563", "GUT_GENOME231247")),]
balloon = ggplot(final_df, aes(x = reorder(Species,Prevalence), y = condition.value, size = Prevalence, fill = log10(Rel_abund))) +
  geom_point(alpha = 0.7, pch=21, colour="black") +
  coord_flip() +
  facet_grid(~ Metadata, scales = "free_x", space = "free") +
  scale_fill_gradient2(low = "navy", mid= "white", high = "red", name=expression(log[10]~(Abundance))) +
  scale_size_continuous(range = c(1,8), name = "Prevalence (%)") +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(panel.spacing = unit(2, "lines")) +
  theme(legend.position = "right") +
  theme(strip.text = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 12, face="italic")) +
  theme(axis.text.x = element_text(angle=45, hjust=0, size = 12))
ggsave(file = "figures/balloon_all.pdf", bg="white", width = 12, height = 4, dpi=300)

# estimate minimum relab
entero.df = abund.data.df[detected.entero,]
entero.relab = abund.data.df/(genome.lengths[rownames(entero.df),"Length"]/1000)
for(col in colnames(entero.relab)) {
  entero.relab[,col] = entero.relab[,col] / sum(entero.relab[,col]) * 100
}
entero.relab[entero.relab == 0] = NA
min(entero.relab, na.rm=TRUE)

# get top E.coli prevalence
ecoli.prev = final_df[which(final_df$Species == "Escherichia coli"),]
ecoli.prev = ecoli.prev[order(ecoli.prev$Prevalence, decreasing=TRUE),]
