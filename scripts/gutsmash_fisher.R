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
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/genofan/")
comb.abund = read.csv("../overlap/combined_coloniz-abund_filt.csv")
rownames(comb.abund) = comb.abund$Species_rep
metadata = comb.abund[,-1]
metadata = metadata[which(metadata$Family != "Enterobacteriaceae"),]
func.data = as.data.frame(fread("gutsmash_results.tsv", header=FALSE))
matrix.df = acast(V1 ~ V3, data=func.data) # V4 class, V3 sublcass

# prepare chisq test
counts.df = merge(matrix.df, metadata, by = "row.names")
counts.df = reshape2::melt(counts.df, id = colnames(counts.df)[which(!colnames(counts.df) %in% colnames(matrix.df))])
counts.df = counts.df[which(counts.df$value !=0),]
counts.df$variable = as.character(counts.df$variable)
counts.df$value = as.numeric(counts.df$value)

# prepare chisq, loop for each variable
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

# clean up labels
barplot.df$variable = gsub("_", " ", barplot.df$variable)
barplot.df$variable = gsub("GR", "Glycyl-radical", barplot.df$variable)
barplot.df$variable = gsub("AA", "amino acid", barplot.df$variable)
barplot.df$variable = gsub("Arginine2", "Arginine to", barplot.df$variable)
barplot.df$variable = gsub("pdu", "Propanediol utilization", barplot.df$variable)
barplot.df$variable = gsub("porA", "R-pyruvate to R-acetate (porA)", barplot.df$variable)
barplot.df$variable = gsub("acetate2butyrate", "Acetate to butyrate", barplot.df$variable)
barplot.df$variable = gsub("succinate2propionate", "Succinate to propionate", barplot.df$variable)
barplot.df$variable = gsub("Others HGD unassigned", "2-hydroxyglutaryl-CoA dehydratase (other)", barplot.df$variable)
barplot.df$variable = gsub("TPP", "Thiamine Pyrophosphate", barplot.df$variable)
barplot.df$variable = gsub("TMA", "Trimethylamine", barplot.df$variable)
barplot.df$variable = gsub("toputrescine", "to putrescine", barplot.df$variable)
barplot.df$variable = gsub("OD", "Oxidative decarboxylation", barplot.df$variable)

bar.plot = ggplot(barplot.df, aes(x = reorder(variable, -value), y = value, fill = type))+
  geom_bar(stat = "identity", position = "stack", alpha=0.8) +
  scale_fill_manual(values = c("tomato", "steelblue"), name="Classification", labels=c("Co-colonizer", "Co-excluder")) +
  coord_flip() +
  theme_classic() +
  ylab("Cluster proportion (%)") +
  theme(legend.text = element_text(size=12), legend.title = element_text(size=12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank())
ggsave(file = "../figures/gutsmash_fisher.pdf", dpi=300, height=5, width=10)
