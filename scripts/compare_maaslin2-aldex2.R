# load libraries
library(plyr)
library(stringr)
library(data.table)
library(matrixStats)
library(ComplexUpset)
library(pheatmap)
library(reshape2)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/overlap")
all.files = list.files()
names(all.files) = all.files
full.dataset = ldply(all.files, read.delim, header=TRUE)

# prepare dataframe
overlap.df = data.frame(matrix(ncol=length(all.files), nrow=length(unique(full.dataset$feature))))
rownames(overlap.df) = unique(full.dataset$feature)
variables = str_split(string = all.files, pattern = "\\.tsv")
colnames(overlap.df) = sapply(variables,"[[",1)

for (col in colnames(overlap.df)){
  subs.df = full.dataset[which(full.dataset$.id == paste(col, ".tsv", sep="")),]
  aldex.col = paste(str_split(subs.df$.id, pattern = "_.*")[[1]][1], "Yes.Est", sep="")
  overlap.df[match(subs.df$feature, rownames(overlap.df)),col] = rowMeans(subs.df[,c("coef", aldex.col)]) # effect size
  #overlap.df[match(subs.df$feature, rownames(overlap.df)),col] = sign(subs.df$coef) # binary
}
overlap.df[is.na(overlap.df)] = 0

# colors for heatmap annotation
annot.df = unique(full.dataset[,c("feature", "Family")])
rownames(annot.df) = annot.df$feature
annot.df = annot.df[,-1, drop=FALSE]
top10.tax = names(sort(table(annot.df[,1]), decreasing=TRUE)[1:10])
annot.df[,1] = ifelse(annot.df[,1] %in% top10.tax, annot.df[,1], "Other")
colors = list(Family=c(brewer.pal(10, "Set3"), "grey"))
names(colors$Family) = c(top10.tax, "Other")

# filter data
overlap.df = overlap.df[,grep("presence", colnames(overlap.df))]

# plot heatmap
pheatmap(t(overlap.df), show_colnames = FALSE, annotation_col = annot.df, annotation_colors = colors)
