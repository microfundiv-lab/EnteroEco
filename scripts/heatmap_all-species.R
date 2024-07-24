# load libraries
library(plyr)
library(stringr)
library(data.table)
library(matrixStats)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ComplexHeatmap)
library(ggpubr)
library(scales)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/overlap")
metadata = read.csv("combined_coloniz-abund_all.csv", check.names = FALSE)
metadata = metadata[order(metadata$ES, decreasing=TRUE),]
rownames(metadata) = metadata$Species_rep

# set annotation and color for Family
annot.df = metadata[,c("Species_rep", "Order")]
rownames(annot.df) = annot.df$Species_rep
annot.df = annot.df[,-1, drop=FALSE]
top.taxa = names(sort(table(metadata$Order), decreasing = TRUE)[1:8])
annot.df$Order = ifelse(annot.df$Order %in% top.taxa, annot.df$Order, "Other")
colors.map = c(brewer.pal(length(top.taxa), "Set3"), "lightgrey")
names(colors.map) = c(top.taxa, "Other")

# color for value, by using 0 as the midpoint
min_val = min(metadata[,c(2:13)], na.rm = TRUE)
max_val = max(metadata[,c(2:13)], na.rm = TRUE)
total_range = max_val - min_val
white_position = abs(min_val) / total_range
colors_before_white = colorRampPalette(c("steelblue", "grey95"))(floor(white_position * 100))
colors_after_white = colorRampPalette(c("grey95", "tomato"))(100 - floor(white_position * 100))
colors_value = c(colors_before_white, colors_after_white)

# rename columns
heatmap.df = t(metadata[,c(2:13)])
rownames(heatmap.df) = gsub("Entero.", "Enterobacteriaceae ", rownames(heatmap.df))
rownames(heatmap.df) = gsub("Ecoli.", "E. coli ", rownames(heatmap.df))
rownames(heatmap.df) = gsub("Kpneumo.", "K. pneumoniae ", rownames(heatmap.df))
rownames(heatmap.df) = gsub("_All", " (All)", rownames(heatmap.df))
rownames(heatmap.df) = gsub("_Healthy-Adult", " (Healthy adults)", rownames(heatmap.df))

# plot heatmap
pdf("../figures/heatmap_all-species.pdf", width=9, height=4)
ha = HeatmapAnnotation(df = annot.df, col = list(Order = colors.map), annotation_name_gp = gpar(fontsize = 10), annotation_legend_param = list(nrow=2))
ht = Heatmap(heatmap.df, show_column_names = FALSE,
              show_row_names = TRUE, top_annotation = ha,
              col = colors_value, cluster_rows = TRUE, heatmap_legend_param = list(title = "Effect size", legend_direction = "horizontal"),
              cluster_columns = FALSE, row_km = 4, row_title = NULL, row_names_gp = gpar(fontsize = 10))
draw(ht, heatmap_legend_side="top", annotation_legend_side="top", padding = unit(c(2, 2, 2, 40), "mm"))
dev.off()
