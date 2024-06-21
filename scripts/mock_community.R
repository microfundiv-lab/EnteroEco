# load libraries
library(ape)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(vegan)
library(ComplexUpset)
library(pheatmap)

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/bwa/mock/")
mock.samples = c("S15", "S16", "S33")
abund.data = data.frame(fread("mock_counts.csv"))
rownames(abund.data) = abund.data$Genome
abund.data = abund.data[,mock.samples]
colnames(abund.data) = c("Mock1", "Mock2", "Mock3")
genomes.metadata = read.delim("../../metadata/genomes_uhgg-v1.2.tsv", stringsAsFactors = FALSE)
rownames(genomes.metadata) = genomes.metadata$Genome
genomes.metadata = genomes.metadata[rownames(abund.data),]

# load taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.df = separate(data = genomes.metadata, col = Lineage, sep = ";", into = ranks)

# calculate relative abundance
bwa.relab = abund.data/(tax.df$Length/1000)
for(col in colnames(bwa.relab)) {
  bwa.relab[,col] = bwa.relab[,col] / sum(bwa.relab[,col]) * 100
}
bwa.relab = bwa.relab[which(rowSums(bwa.relab > 0) > 0),]
bwa.relab$Expected = 100/8

# prepare df for plot
melt.df = reshape2::melt(as.matrix(bwa.relab))
melt.df = melt.df[which(melt.df$value > 0),]
colnames(melt.df)[1] = "Genome"
melt.df = merge(melt.df, tax.df, by="Genome")
melt.df$Species = gsub("s__", "", melt.df$Species)

# stacked bar mock community
mock.plot = ggplot(melt.df, aes(x=Var2, y=value, fill=Species)) +
  geom_bar(stat="identity") +
  theme_classic() +
  xlab("") +
  ylab("Relative abundance (%)") +
  scale_x_discrete(limits=c("Expected", "Mock1", "Mock2", "Mock3")) +
  theme(axis.text = element_text(size=14)) +
  theme(axis.title = element_text(size=14)) +
  theme(legend.text = element_text(size=12, face="italic")) +
  theme(legend.title = element_text(size=12)) +
  scale_fill_manual(values=c(brewer.pal(12, "Set3")))
ggsave(mock.plot, filename="../../figures/stacked_mock.pdf", width=7, height=5)
