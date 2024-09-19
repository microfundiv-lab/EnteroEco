# load libraries
library(data.table)
library(pheatmap)
library(tidyr)

# load input
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
comb.abund = read.csv("overlap/combined_coloniz-abund_all.csv")
abund.data = as.data.frame(fread("bwa/bwa_counts-filtered_samples.csv"))
rownames(abund.data) = abund.data$Genome
abund.data = abund.data[,-1]

# parse taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
species.data = read.delim("metadata/species_uhgg_v1.2.tsv")
species.data = unique(species.data[,c("Species_rep", "Lineage")])
species.data = separate(data = species.data, col = Lineage, sep = ";", into = ranks)
rownames(species.data) = species.data$Species_rep
species.data = data.frame(lapply(species.data, function(x) { gsub(".__", "", x)}))
rownames(species.data) = species.data$Species_rep

# select top taxa
n = 1000
top.taxa = names(sort(rowSums(abund.data > 0), decreasing=TRUE)[1:n])
species.data.top = species.data[top.taxa,]
species.data.top$Classification = "Not significant"
co.excluders = comb.abund[which(comb.abund$Classification == "Co-excluder"), "Species_rep"]
co.colonizers = comb.abund[which(comb.abund$Classification == "Co-colonizer"), "Species_rep"]
species.data.top[intersect(co.excluders, top.taxa), "Classification"] = "Co-excluder"
species.data.top[intersect(co.colonizers, top.taxa), "Classification"] = "Co-colonizer"

# plot
order.props = as.data.frame.matrix(table(species.data.top[,c("Order", "Classification")]))
order.props.rows = cbind(order.props/rowSums(order.props))*100
order.props.cols = 100 * order.props/colSums(order.props)[col(order.props)]
pheatmap(log10(order.props+1), angle_col=45, legend_labels=c("Test"),
         filename = "figures/heatmap_top-prev.pdf", width=4, height=9)
under.taxa = order.props[which(order.props$`Not significant` > 0 & order.props$`Co-colonizer` == 0 & order.props$`Co-excluder` == 0),]
under.taxa = under.taxa[order(under.taxa$`Not significant`, decreasing=TRUE),]
