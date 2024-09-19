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
metadata = metadata[which(metadata$Phylum %in% c("Bacillota", "Bacillota_A", "Bacillota_C")),]

# parse COG results
cog.data = fread("cog_results.tsv")
cog.data = cog.data[which(cog.data$COG_category != "-"),]
cog.data = separate_rows(cog.data, COG_category, sep="(?!^)")
cog.data = as.data.frame(cog.data[which(cog.data$COG_category != ""),])
cog.data$`#query` = sub("_\\d+$", "", cog.data$`#query`)
cog.matrix = acast(COG_category ~ `#query`, data=cog.data)
cog.catalog = openxlsx::read.xlsx('COG_categories.xlsx')
rownames(cog.catalog) = cog.catalog$variable
cog.relab = t(cog.matrix)/colSums(cog.matrix)*100
cog.meta = merge(cog.relab, metadata, by="row.names")

# wilcox test
cogs = colnames(cog.relab)
wilcox.res = lapply(cogs, function(x) {
  mean.pos = mean(cog.meta[cog.meta$Classification == "Co-colonizer", x])
  mean.neg = mean(cog.meta[cog.meta$Classification == "Co-excluder", x])
  wilcox.formula = formula(paste(paste0(x, " ~ Classification")))
  effsize = wilcox_effsize(cog.meta, wilcox.formula)$effsize
  effsize = ifelse(mean.pos > mean.neg, effsize, -effsize)
  sign = wilcox_test(cog.meta, wilcox.formula)$p
  return(c(effsize, sign))
})

# process output
raw.output = data.frame(t(data.frame(wilcox.res)))
rownames(raw.output) = cogs
raw.output$FDR = p.adjust(raw.output[,2])
sign.output = raw.output[which(raw.output$FDR < 0.05),]
sign.output = merge(sign.output, cog.catalog, by="row.names")
colnames(sign.output)[1:4] = c("COG", "EffectSize", "Pvalue", "FDR")

# plot barplot
bar.plot = ggplot(sign.output, aes(x = reorder(Categories, EffectSize), y = EffectSize, fill = Group)) +
  geom_bar(stat = "identity", alpha=0.8) +
  scale_fill_manual(values = c("Cellular processes and signaling" = "#C0392B", 
                               "Information storage and processing" = "#56B4E9", 
                               "Metabolism" = "#009E73",
                               "Poorly characterized" = "#979A9A"), name="Functional category")+
  coord_flip() +
  theme_classic() +
  ylab("Effect size") +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "top")
ggsave(bar.plot, file = "../figures/cog_wilcox_bacillota.pdf", dpi=300, height=5, width=7)
