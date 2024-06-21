# load libraries
library(data.table)
library(tidyr)
library(ggplot2)

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/")
abund.data = fread("data/bwa/bwa_counts-filtered_batch-corr_samples.csv")
abund.data.df = data.frame(abund.data[,-1], check.names = FALSE)
rownames(abund.data.df) = abund.data$Genome

# filter metadata
metadata = read.delim("data/metadata/metagenomes_11-2023_samples.tsv", check.names = TRUE)
rownames(metadata) = metadata$Sample
metadata$Disease.name = ifelse(metadata$Health.state == "Healthy", "Healthy", metadata$Disease.name)
abund.data.df = abund.data.df[,metadata$Sample]

# load taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
genomes.metadata = read.delim("data/metadata/species_uhgg_v1.2.tsv")
rownames(genomes.metadata) = genomes.metadata$Genome
tax.df = unique(genomes.metadata[,c("Species_rep", "Lineage")])
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
rownames(tax.df) = tax.df$Species_rep

# define entero presence/absence
entero.species = as.vector(tax.df[which(tax.df$Family == "f__Enterobacteriaceae"),"Species_rep"])
entero.pos = names(which(colSums(abund.data.df[entero.species,]) > 0))
metadata$Entero.presence = ifelse(metadata$Sample %in% entero.pos, "Yes", "No")

# check metadata proportion
entero.absence = metadata[which(metadata$Entero.presence == "No"),]
entero.absence.cont = data.frame(table(entero.absence$Continent))
entero.absence.cont$Metadata = "Continent"
entero.absence.disease = data.frame(table(entero.absence$Disease.name))
entero.absence.disease$Metadata = "Health state"
entero.absence.age = data.frame(table(entero.absence$Age.group))
entero.absence.age$Metadata = "Age group"
entero.absence.df = rbind(entero.absence.cont, entero.absence.disease, entero.absence.age)
entero.absence.df$Entero = "Absence"

entero.presence = metadata[which(metadata$Entero.presence == "Yes"),]
entero.presence.cont = data.frame(table(entero.presence$Continent))
entero.presence.cont$Metadata = "Continent"
entero.presence.disease = data.frame(table(entero.presence$Disease.name))
entero.presence.disease$Metadata = "Health state"
entero.presence.age = data.frame(table(entero.presence$Age.group))
entero.presence.age$Metadata = "Age group"
entero.presence.df = rbind(entero.presence.cont, entero.presence.disease, entero.presence.age)
entero.presence.df$Entero = "Presence"

# plot metadata
entero.df = rbind(entero.absence.df, entero.presence.df)
entero.df$Var1 = gsub("CRC", "Colorectal cancer", entero.df$Var1)
entero.df$Var1 = gsub("CD", "Crohn's disease", entero.df$Var1)
entero.df$Var1 = gsub("UC", "Ulcerative colitis", entero.df$Var1)
entero.df$Var1 = gsub("T1D", "Type 1 diabetes", entero.df$Var1)
entero.df$Var1 = gsub("T2D", "Type 2 diabetes", entero.df$Var1)
entero.df$Var1 = gsub("T2D", "Type 2 diabetes", entero.df$Var1)
entero.df$Var1 = gsub("AS", "Ankylosing spondylitis", entero.df$Var1)
entero.df$Var1 = gsub("RA", "Rheumatoid arthritis", entero.df$Var1)
entero.df$Var1 = gsub("Parkinson", "Parkinson's disease", entero.df$Var1)
entero.df$Var1 = gsub("BoneDisease", "Bone disease", entero.df$Var1)
entero.df$Var1 = gsub("MS", "Multiple sclerosis", entero.df$Var1)
entero.df$Var1 = gsub("ME", "Myalgic encephalomyelitis", entero.df$Var1)
entero.df$Var1 = factor(entero.df$Var1, levels=sort(unique(entero.df$Var1), decreasing=TRUE))

metadata.plot = ggplot(entero.df, aes(x=Var1, y=log10(Freq), fill=Entero)) +
  geom_bar(stat="identity", position="dodge", alpha=0.8, width=0.5) +
  facet_wrap(~ Metadata, scales = "free_y") +
  scale_fill_manual(values=c("steelblue", "tomato"), name="Colonization status") +
  theme_classic() +
  theme(legend.title = element_text(size=16)) +
  theme(legend.text = element_text(size=14)) +
  theme(strip.background = element_blank()) +
  theme(strip.text = element_text(face = "bold", size=16)) +
  theme(axis.title = element_text(size=16)) +
  theme(axis.text = element_text(size=14)) +
  theme(legend.position = "bottom") +
  coord_flip() +
  xlab("") +
  ylab(expression(log[10]~(Number~of~Samples)))
ggsave(file="data/figures/metadata_distr.pdf", height=5, width=10)

# calculate proportions
entero.presence.df$Prop = entero.presence.df$Freq/(entero.presence.df$Freq+entero.absence.df$Freq)*100
median.cont = quantile(entero.presence.df$Prop[which(entero.presence.df$Metadata == "Continent")])
median.age = quantile(entero.presence.df$Prop[which(entero.presence.df$Metadata == "Age group")])
median.health = quantile(entero.presence.df$Prop[which(entero.presence.df$Metadata == "Health state")])
