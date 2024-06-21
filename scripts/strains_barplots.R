# load libraries
library(reshape2)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(rstatix)
library(ggpubr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/strains/")
all.cog = fread("all_genes_species_cogs.tsv", header = F)
all.cog$Group = "All"
sig.cog = fread("sig_genes_species_cogs.tsv", header = F)
sig.cog$Group = "Significant"
cog.data = rbind(all.cog, sig.cog)
cog.data = cog.data[,-1]
colnames(cog.data)[1:2] = c("#query", "COG_category")
cog.data = as.data.frame(cog.data[which(cog.data$COG_category != ""),])
cog.data$value = 1

# prepare chisq, loop for each variable
result.df = data.frame(variable=character(), sig = numeric(), all = numeric(), pvalue=numeric(), result=numeric())
var.list = unique(cog.data$COG_category)
for (variable in var.list){
  cog.data$variable.group = ifelse(cog.data$COG_category == variable, variable, "all_other")
  cont.table = with(cog.data, xtabs(value ~ Group + variable.group)) # create contingency table
  test.result = fisher.test(cont.table)
  sig_prop = cont.table["Significant",variable]/(cont.table["Significant","all_other"]+cont.table["Significant",variable])*100
  all_prop = cont.table["All",variable]/(cont.table["All","all_other"]+cont.table["All",variable])*100
  result = sig_prop - all_prop
  result.df = rbind(result.df, data.frame(variable = variable, sig = sig_prop, all = all_prop, pvalue = test.result$p.value, result = result))
  cog.data$variable.group = NULL # clear for next cycle
}
result.df$FDR = p.adjust(result.df$pvalue, method="fdr")
result.df = result.df[which(result.df$FDR < 0.05),]
result.df$Classification = ifelse(result.df$result > 0, "Higher", "Lower")
rownames(result.df) = NULL

# prepare data for bar plot
barplot.df = reshape2::melt(result.df)
colnames(barplot.df) = c("variable", "Classification", "type", "value")
barplot.df = barplot.df[which(barplot.df$type %in% c("sig", "all")),]

# add category annotation
cog.catalog = openxlsx::read.xlsx('../genofan/COG_categories.xlsx')
rownames(cog.catalog) = cog.catalog$variable
order.cats = unique(barplot.df[order(barplot.df$value, decreasing=TRUE),"variable"])
order.annot = cog.catalog[order.cats,"Categories"]
barplot.df$Cat = cog.catalog[match(barplot.df$variable, cog.catalog$variable),"Group"]

# plot barplot
bar.plot = ggplot(barplot.df, aes(x = variable, y = value, fill = type))+
  geom_bar(stat = "identity", position = "dodge", width=0.8, alpha=0.6) +
  scale_fill_manual(values = c("purple4", "darkgrey"), name="Classification", labels=c("Significant", "All")) +
  scale_x_discrete(limits=order.cats, labels=c("Function\nunknown", "Replication\nand repair", "Translation", "Nucleotide\nmetabolism")) +
  theme_classic() +
  ylab("Cluster proportion (%)") +
  theme(legend.position = "right", legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_blank())

# prepare gene count data
gene.counts = read.delim("genes_per_species.txt", sep=" ", header=FALSE)
rownames(gene.counts) = gsub("/glm-filter.log:Found", "", gene.counts$V1)
gene.counts = gene.counts[,"V2", drop=FALSE]
gene.counts = gene.counts[order(gene.counts[,1], decreasing=TRUE),,drop=FALSE]

# parse taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
species.metadata = read.delim("../metadata/species_uhgg_v1.2.tsv", stringsAsFactors = FALSE)
tax.df = unique(species.metadata[,c("Species_rep", "Lineage", "Status")])
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
tax.df$Taxon = NA
for (n in 1:nrow(tax.df)){
  for (r in ranks){
    if(!grepl("__$", tax.df[n,r])){
      tax.df[n,"Taxon"] = tax.df[n,r] }}}
tax.df = data.frame(lapply(tax.df, function(x) { gsub(".__", "", x)}))
rownames(tax.df) = tax.df$Species_rep
gene.counts$Taxon = tax.df[rownames(gene.counts), "Taxon"]
gene.counts = rbind(gene.counts, c(sum(gene.counts[3:nrow(gene.counts),"V2"]), "Other"))
rownames(gene.counts)[nrow(gene.counts)] = "Other"
gene.counts$V2 = as.numeric(gene.counts$V2)

# gene counts per cog
selected.taxa = c("GUT_GENOME096210", "GUT_GENOME000718")
#sig.cog$V2 = ifelse(sig.cog$V2 %in% selected.taxa, sig.cog$V2, "Other")
sig.cog$Cat = cog.catalog[match(sig.cog$V3, cog.catalog$variable),"Group"]
cog.counts = aggregate(V1 ~ V2 + Cat, data=sig.cog, FUN=length)
for (species in unique(cog.counts$V2)) {
  n_annot = sum(cog.counts[which(cog.counts$V2 == species),"V1"])
  n_none = gene.counts$V2[which(rownames(gene.counts) == species)] - n_annot
  cog.counts = rbind(cog.counts, c(species, "Unassigned", n_none))
  cog.counts$V1 = as.numeric(cog.counts$V1)
}
cog.counts$Taxon = gene.counts[match(cog.counts$V2, rownames(gene.counts)),"Taxon"]

# plot counts
gene.plot = ggplot(cog.counts, aes(x=reorder(Taxon, V1), y=V1, fill=Cat)) +
  geom_bar(stat="identity", alpha=0.6, width=0.8) +
  theme_classic() +
  coord_flip() +
  xlab("") +
  ylab("# Significant accessory genes") +
  scale_fill_manual(values=c("tomato", "steelblue", "darkgreen", "darkgrey", "lightgrey"), name="Functional category") +
  theme(legend.position = "inside", legend.position.inside = c(0.6,0.45), legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12, face="italic"))

# combine plots
ggarrange(gene.plot, bar.plot, widths=c(1.2,1), labels=c("a", "b"), font.label = list(size=18))
ggsave(file = "../figures/strains_barplots.pdf", dpi=300, height=5, width=14)
