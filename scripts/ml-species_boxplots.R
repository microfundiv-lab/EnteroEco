# load libraries
library(data.table)
library(tidyverse)
library(ggplot2)

# setwd
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/ml-species/")

# load performance data
input.files =  list.files(pattern = "results.csv", full.names=TRUE)
pf.list = lapply(input.files, function(i) {
  analysis = strsplit(i, "/")[[1]][2]
  taxon = strsplit(analysis, "_")[[1]][2]
  meta = strsplit(analysis, "_")[[1]][1]
  df = read.csv(i); df$Analysis = analysis; df$Taxon = taxon; df$Dataset = meta
  return(df)
})
pf.combined = as.data.frame(rbindlist(pf.list))

# select best
pf.agg = aggregate(AUC ~ method + Taxon, data=pf.combined, FUN=median)
pf.agg = pf.agg[order(pf.agg$AUC, decreasing=TRUE),]

# pf boxplots
pf.combined$method = gsub("glmnet", "Ridge Regression", pf.combined$method)
pf.combined$method = gsub("rf", "Random Forest", pf.combined$method)
pf.combined$method = gsub("xgbTree", "Gradient Boosting", pf.combined$method)
pf.combined$method = factor(pf.combined$method, levels=c("Gradient Boosting", "Random Forest", "Ridge Regression"))
pf.combined$Taxon = factor(pf.combined$Taxon, levels=c("Entero", "Ecoli", "Kpneumo"))
ml.compare = ggplot(pf.combined, aes(x=Dataset, y=AUC, fill=Taxon)) +
  geom_point(alpha=0.6, size=0.9, colour="darkgrey", position = position_jitterdodge(jitter.width=0.2, dodge.width = 0.7)) +
  geom_boxplot(alpha=0.5, outlier.colour = NA, colour="black", width=0.7) +
  facet_wrap(~ method, ncol=3) +
  geom_hline(yintercept = 0.7, linetype="dashed") +
  ylim(0.5,1) +
  ylab("AUC") +
  theme_bw() +
  scale_x_discrete(limits=c("HQ", "Healthy-Adult"), labels=c("All", "Healthy adults")) +
  scale_fill_manual(values=c("steelblue", "darkgreen", "tomato"), labels=c("Enterobacteriaceae", "E. coli", "K. pneumoniae")) +
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.pos = "top") +
  theme(strip.background = element_blank(), strip.text = element_text(size=14)) +
  theme(legend.position="top", legend.box = "horizontal", legend.text=element_text(size=14, face="italic"),
        legend.title = element_blank()) +
  theme(panel.grid.minor = element_blank()) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.y = element_text(size=14)) + 
  theme(axis.title.y = element_text(size=16)) + 
  theme(axis.text.x = element_text(size=14))
ggsave(filename="../figures/ml-species_boxplots.pdf", width=8, height=4)
