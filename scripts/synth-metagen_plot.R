# load libraries
library(ggplot2)
library(tidyr)
library(grid)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/synth_metagen/")
metadata = read.delim("../metadata/metagenomes_synth.tsv")
abund.data = read.csv("../bwa/synth_metagen/bwa_counts-filtered.csv")
rownames(abund.data) = abund.data$Genome
abund.data.df = abund.data[,-1]
genomes.metadata = read.delim("../metadata/species_uhgg_v1.2.tsv")

# load taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(genomes.metadata) = genomes.metadata$Genome
tax.df = separate(data = genomes.metadata, col = Lineage, sep = ";", into = ranks)

# calculate relative abundance
abund.data.df = abund.data.df[tax.df$Species_rep,]
bwa.relab = abund.data.df/(tax.df$Length/1000)
for(col in colnames(bwa.relab)) {
  bwa.relab[,col] = bwa.relab[,col] / sum(bwa.relab[,col]) * 100
}

# reshape
bwa.melt = reshape2::melt(cbind(rownames(bwa.relab), bwa.relab))
bwa.melt$Sample = as.vector(bwa.melt$variable)
bwa.fi = merge(bwa.melt, metadata, by="Sample")
bwa.fi = bwa.fi[which(bwa.fi$`rownames(bwa.relab)` == bwa.fi$Species),]
bwa.fi$Depth = factor(bwa.fi$Depth)

#
count_above_zero = function(x) {sum(x > 0)}

# aggregate to get mean and sd
bwa.mean = aggregate(value ~ Depth + Abundance, data=bwa.fi, FUN=mean)
bwa.sd = aggregate(value ~ Depth + Abundance, data=bwa.fi, FUN=sd)
bwa.count = aggregate(value ~ Depth + Abundance, data=bwa.fi, FUN=function(x) {sum(x > 0)})
bwa.mean$sd = bwa.sd$value
bwa.mean$Detected = bwa.count$value
bwa.mean$Depth_class = ifelse(bwa.mean$Depth == 12640000, "Low depth (13M)", ifelse(bwa.mean$Depth == 30869000, "Medium depth (31M)", "High depth (51M)"))
bwa.mean$Depth_class = factor(bwa.mean$Depth_class, levels=c("Low depth (13M)", "Medium depth (31M)", "High depth (51M)"))
bwa.mean$LOD = ifelse(bwa.mean$Depth_class == "Low depth (13M)", 0.01, ifelse(bwa.mean$Depth_class == "Medium depth (31M)", 0.005, 0.003))
bwa.mean = bwa.mean[which(bwa.mean$Detected == 5),]

# plot
lod.plot = ggplot(bwa.mean, aes(x=Abundance, y=value)) +
  geom_point(size=1.5, colour="black") +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.1) +
  geom_hline(data = bwa.mean, aes(yintercept=LOD), linetype="dashed") +
  facet_wrap(~ Depth_class) +
  theme_bw() +
  xlab(bquote("Spiked-in abundance (%)")) +
  ylab(bquote("Detected abundance (%)")) +
  scale_x_log10(limits=c(0.001,1.5), breaks=c(0.001,0.003,0.005,0.01,0.05,0.1,0.5,1.0), labels=c(0.001,0.003,0.005,0.01,0.05,0.1,0.5,1.0)) +
  scale_y_log10(limits=c(0.001,1.5), breaks=c(0.001,0.003,0.005,0.01,0.05,0.1,0.5,1.0), labels=c(0.001,0.003,0.005,0.01,0.05,0.1,0.5,1.0)) +
  guides(colour="none") +
  theme(strip.background = element_blank()) +
  theme(strip.text = element_text(size=14)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1))
ggsave(filename = "../figures/synth-metagen_plot.pdf", dpi=300, width=12, height=4)
