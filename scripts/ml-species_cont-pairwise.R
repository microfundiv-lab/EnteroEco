# load libraries
library(tidyverse)
library(data.table)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)

# set working directory
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/ml-species/continent_pairwise/")

# load performance data
input.files =  list.files(pattern = "results.csv", full.names=TRUE)
pf.list = lapply(input.files, function(i) {
  filename = strsplit(i, "/")[[1]][2]
  continent1 = strsplit(filename, "_")[[1]][1]
  continent2 = strsplit(filename, "_")[[1]][3]
  df = read.csv(i); df$Training = continent1; df$Testing = continent2
  return(df)
})

# combine all results into a dataframe
pf.combined = as.data.frame(rbindlist(pf.list))
pf.agg = aggregate(AUC ~ Training + Testing, FUN=median, data=pf.combined)

# rename diseases
rename_disease = function(df, column) {
  df[,column] = gsub("NorthAmerica", "North America", df[,column])
  df[,column] = gsub("SouthAmerica", "South America", df[,column])
  df2 = df
  return(df2)
}

pf.agg.ren1 = rename_disease(pf.agg, "Training")
pf.agg.fi = rename_disease(pf.agg.ren1, "Testing")
pf.agg.fi$AUC = round(pf.agg.fi$AUC, digits=2)
pf.agg.fi$Class = ifelse(pf.agg.fi$AUC >= 0.8, ">=0.8", ifelse(pf.agg.fi$AUC >= 0.7, "0.7-0.79", ifelse(pf.agg.fi$AUC > 0.6, "0.6-0.69", ifelse(pf.agg.fi$AUC >0.5, "0.5-0.59", "<0.5"))))
pf.agg.fi$Class = factor(pf.agg.fi$Class, levels = c("<0.5", "0.5-0.59", "0.6-0.69", "0.7-0.79", ">=0.8"))


# plot heatmap
ml.heat = ggplot(pf.agg.fi, aes(x=Training, y=Testing, fill=Class, label=AUC)) +
  geom_tile() +
  geom_text() +
  theme_classic() +
  scale_fill_manual(values=c("grey90", "lightblue", "steelblue", "pink2", "tomato"), name="AUROC") +
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14))
ggsave(file="../../figures/ml_cont-pairwise.pdf", dpi=300, width=6, height=4)
