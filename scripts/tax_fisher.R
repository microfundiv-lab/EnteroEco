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
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data")
comb.abund = read.csv("overlap/combined_coloniz-abund_filt.csv")
comb.abund = comb.abund[which(comb.abund$Family != "Enterobacteriaceae"),]
rownames(comb.abund) = comb.abund$Species_rep
metadata = comb.abund[,-1]

# select rank
rank = "Order"
tax.df = metadata[,rank, drop=FALSE]
tax.df$Genome = rownames(tax.df)
tax.df$Count = 1
matrix.df = acast(as.formula(paste("Genome", rank, sep="~")), value.var = "Count", data=tax.df)
matrix.df[is.na(matrix.df)] = 0

# check correlation
cor.df = as.data.frame.matrix(table(comb.abund[,c(rank, "Classification")]))
cor = cor.test(cor.df[,1], cor.df[,2])

# prepare fisher test
counts.df = merge(matrix.df, metadata, by = "row.names")
counts.df = reshape2::melt(counts.df, id = colnames(counts.df)[which(!colnames(counts.df) %in% colnames(matrix.df))])
counts.df = counts.df[which(counts.df$value !=0),]
counts.df$variable = as.character(counts.df$variable)
counts.df$value = as.numeric(counts.df$value)

# prepare fisher, loop for each variable
result.df = data.frame(variable=character(), pos = numeric(), neg = numeric(), pvalue=numeric(), result=numeric())
var.list = unique(counts.df$variable)
for (variable in var.list){
  counts.df$variable.group = ifelse(counts.df$variable == variable, variable, "all_other")
  cont.table = with(counts.df, xtabs(value ~ Classification + variable.group)) # create contingency table
  test.result = fisher.test(cont.table)
  pos_prop = cont.table["Co-colonizer",variable]/(cont.table["Co-colonizer","all_other"]+cont.table["Co-colonizer",variable])*100
  neg_prop = cont.table["Co-excluder",variable]/(cont.table["Co-excluder","all_other"]+cont.table["Co-excluder",variable])*100
  result = pos_prop - neg_prop
  result.df = rbind(result.df, data.frame(variable = variable, pos = pos_prop, neg = neg_prop, pvalue = test.result$p.value, result = result))
  counts.df$variable.group = NULL # clear for next cycle
}
result.df$FDR = p.adjust(result.df$pvalue, method="fdr")
result.df = result.df[which(result.df$FDR < 0.05),]
result.df$Classification = ifelse(result.df$result > 0, "Co-colonizer", "Co-excluder")
rownames(result.df) = NULL

# bar plots
barplot.df = reshape2::melt(result.df)
colnames(barplot.df) = c("variable", "Classification", "type", "value")
barplot.df = barplot.df[which(barplot.df$type %in% c("pos", "neg")),]
bar.plot = ggplot(barplot.df, aes(x = reorder(variable, -value), y = value, fill = type))+
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("tomato", "steelblue"), name="Classification", labels=c("Co-colonizer", "Co-excluder")) +
  coord_flip() +
  theme_classic() +
  ylab("Taxon proportion (%)") +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank())
ggsave(file = "figures/tax_fisher.png", dpi=300, height=4, width=7)
