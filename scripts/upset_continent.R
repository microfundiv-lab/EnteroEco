# load libraries
library(data.table)
library(matrixStats)
library(ComplexUpset)
library(reshape2)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")

#load species info, match to the enterobacteriaceae to sum
genomes.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
rownames(genomes.metadata) = genomes.metadata$Genome
genome.lengths = genomes.metadata[unique(as.vector(genomes.metadata$Species_rep)),c("Genome", "Length")]
gtdb.tax = unique(genomes.metadata[,c("Species_rep", "Lineage", "Genome_type", "Status")])
rownames(gtdb.tax) = gtdb.tax$Species_rep
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
gtdb.tax = separate(data = gtdb.tax, col = Lineage, sep = ";", into = ranks)
colnames(gtdb.tax)[1]= "Genome"
gtdb.tax = gtdb.tax[,-1]

# load data
abund.data = fread("bwa/bwa_counts-filtered_batch-corr_samples.csv") 
abund.data.df = data.frame(abund.data[,-1], check.names = FALSE)
rownames(abund.data.df) = abund.data$Genome
da.result = abund.data.df
da.entero = da.result[rownames(da.result) %in% rownames(gtdb.tax[gtdb.tax$Family == "f__Enterobacteriaceae",]),]
da.entero = da.entero[rowSums(da.entero) > 0, ]

#filter metadata for different variables 
metadata = read.delim("metadata/metagenomes_11-2023_samples.tsv", stringsAsFactors = FALSE, check.names = TRUE) #, na.strings = ""
rownames(metadata) = metadata$Sample
metadata = metadata[which(!is.na(metadata$Disease.name) & !is.na(metadata$Age.group) & !is.na(metadata$Continent)),]

# change metadata here
meta.filt = metadata[,c("Sample", "Read.count","Age.group", "Disease.name", "Country", "Continent")]
rownames(meta.filt) = meta.filt$Sample

# load Entero species info
entero.list = readLines("metadata/species_entero.txt")

# select species, just select 5 species
abund.data.entero = abund.data.df[entero.list,]
ranked.prev = data.frame(sort(rowSums(abund.data.entero > 0), decreasing=TRUE))
ranked.prev.tax = merge(ranked.prev, gtdb.tax, by="row.names")
colnames(ranked.prev.tax)[2] = "Prevalence"
ranked.prev.tax = ranked.prev.tax[order(ranked.prev.tax$Prevalence, decreasing=TRUE),c("Row.names", "Prevalence", "Species")]

da.sign = da.entero[rownames(da.entero) %in% entero.list,]
da.sign = da.entero[rownames(da.entero) %in% c("GUT_GENOME144544", "GUT_GENOME147598","GUT_GENOME143760", "GUT_GENOME000563", "GUT_GENOME231247"),]
da.sign[da.sign != 0] = 1
da.sign = as.data.frame(t(da.sign))

# rename colnames, top5 genus from ranked.prev.tax
colnames(da.sign) = gsub("GUT_GENOME144544", "Escherichia coli", colnames(da.sign))
colnames(da.sign) = gsub("GUT_GENOME147598", "Klebsiella pneumoniae", colnames(da.sign))
colnames(da.sign) = gsub("GUT_GENOME143760", "Enterobacter hormaechei_A", colnames(da.sign))
colnames(da.sign) = gsub("GUT_GENOME000563", "Citrobacter freundii", colnames(da.sign))
colnames(da.sign) = gsub("GUT_GENOME231247", "Kluyvera ascorbata", colnames(da.sign))

# meger metadate with sample info
da.select = merge(da.sign, meta.filt, by = "row.names")

# generate upset
variable = "Continent"
upset.df = da.select
upset.df = upset.df[rowSums(upset.df[,2:6]) > 0,] #for species
keep.taxa = names(sort(table(upset.df[,variable]), decreasing=TRUE)[1:10])
upset.df[,variable] = ifelse(upset.df[,variable] %in% keep.taxa, upset.df[,variable], "Other")
upset.df[,variable] = factor(upset.df[,variable], levels=c(sort(keep.taxa), "Other"))

upset.plot = upset(upset.df, colnames(upset.df)[2:6], n_intersections=15, name="",
                   base_annotations=list('Intersection size'=intersection_size(
                     mapping=aes(fill=upset.df[[variable]]))
                     + scale_fill_manual(values=c(brewer.pal(10, "Set3"), "darkgrey"), name=variable)),
                   annotations = list(
                     variable=(
                       ggplot(mapping=aes(fill=upset.df[[variable]])) 
                       + geom_bar(stat='count', position='fill')
                       + scale_y_continuous(labels=scales::percent_format())
                       + scale_fill_manual(values=c(brewer.pal(10, "Set3"), "darkgrey"), name=variable)
                       + ylab('% of samples'))),
                   set_sizes = upset_set_size() + ylab('Number of samples'),
                   width_ratio=0.25,
                   themes=upset_default_themes(text=element_text(size=14)))
ggsave(upset.plot, filename="figures/upset_continent.pdf", height=8, width=12, dpi=300)

