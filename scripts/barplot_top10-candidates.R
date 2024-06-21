# load libraries
library(plyr)
library(stringr)
library(data.table)
library(matrixStats)
library(ComplexUpset)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(ggpubr)
library(scales)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/overlap")
metadata = read.csv("combined_coloniz-abund_filt.csv", check.names = FALSE)

top.negative = metadata$Species_rep[order(metadata$`Co-excluder`, -metadata$ES, decreasing=TRUE)][1:10]
top.positive = metadata$Species_rep[order(metadata$`Co-colonizer`, metadata$ES, decreasing=TRUE)][1:10]
top.df = metadata[which(metadata$Species_rep %in% top.positive | metadata$Species_rep %in% top.negative),]

# rename species if unassigned
top.df$Species = ifelse(is.na(top.df$Species) | top.df$Species == "", paste(top.df$Genus, "sp."), top.df$Species)

# set colors
colors = brewer.pal(length(unique(top.df$Family)), "Set1")
names(colors) = unique(top.df$Family)

# selected columns
data.vars = c("Ecoli.colonization_All", "Ecoli.abundance_All", "Kpneumo.colonization_All", 
              "Kpneumo.abundance_All", "Entero.colonization_All", "Entero.abundance_All",
              "Ecoli.colonization_Healthy-Adult", "Ecoli.abundance_Healthy-Adult", "Kpneumo.colonization_Healthy-Adult", 
              "Kpneumo.abundance_Healthy-Adult", "Entero.colonization_Healthy-Adult", "Entero.abundance_Healthy-Adult")

# bar plot
top.df[top.df == 0] = NA
top.barplot = top.df %>%
  mutate(sd = apply(.[,data.vars], 1, sd, na.rm = TRUE)) %>% 
  ggplot(aes(y = ES, x = reorder(Species, ES), fill = Family)) +
  geom_col(position = "dodge", alpha = 0.5)+
  scale_fill_manual(values = colors) +
  coord_flip() + 
  theme_classic() +
  ylab("Effect size") +
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size = 12, face = 'italic'), axis.title.y = element_blank(), axis.title.x = element_text(size=14),
        legend.text = element_text(size = 12, face = 'italic', hjust = 0), legend.title = element_text(size=12), legend.position = "right")
ggsave(top.barplot, filename = "../figures/barplot_top10-candidates.pdf", height = 5, width = 7, dpi = 300)