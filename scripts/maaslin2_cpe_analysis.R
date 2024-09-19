# load libraries
library(Maaslin2)
library(tidyr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
metadata = read.delim("metadata/metagenomes_timeseries.tsv")
abund.data = read.csv("bwa/timeseries/bwa_counts-filtered.csv")
genomes.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
output_dir = "timeseries/maaslin2"

# format data
rownames(metadata) = metadata$Run
abund.data.df = data.frame(abund.data[,-1], check.names = FALSE)
rownames(abund.data.df) = abund.data$Genome
abund.data.df = abund.data.df[,metadata$Run]
abund.data.df = abund.data.df[,which(colSums(abund.data.df) > 0)]
metadata = metadata[colnames(abund.data.df),]
metadata = metadata[which(metadata$Antibiotics.since.last.visit != "Yes" & metadata$Hospitalization.since.last.visit != "Yes"),]
metadata$CPE_status = ifelse(metadata$Colonization.status...group == "CPE_positive", "CPE_positive", 
                             ifelse(metadata$Colonization.status...group %in% c("CPE_negative_0", "CPE_negative_1", "CPE_negative_2"), 
                                   "CPE_negative_index", "CPE_negative_control"))

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

# run Maaslin2
fit_data = Maaslin2(cores = 4, normalization = "NONE", transform = "LOG",
  min_prevalence = 0.01, max_significance = 0.05, correction = "BH",
  input_data = bwa.relab, input_metadata = metadata, output = output_dir,
  fixed_effects = c("CPE_status"),
  reference = c("CPE_status,CPE_positive"),
  random_effects = "Individual.Code",
  plot_heatmap = FALSE, plot_scatter = FALSE)