# load libraries
library(RColorBrewer)
library(maptools)
library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(grid)
library(rgeos)
library(maptools)
library(mapproj)

# load input
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
metadata = read.delim("metadata/metagenomes_07-2024_samples.tsv")
metadata = metadata[!is.na(metadata$Country),]
rownames(metadata) = metadata$Sample
metadata$Country = gsub("United States of America", "United States", metadata$Country)

# filter metadata
metadata = metadata[which(!is.na(metadata$Disease.name) & !is.na(metadata$Age.group) & !is.na(metadata$Continent)),]

# count samples per country
countries = unique(metadata$Country)

ddf = data.frame(matrix(NA, ncol=2, nrow=length(countries)), row.names=countries)
colnames(ddf) = c("country", "samples")
ddf$country = rownames(ddf)

# count metagenomes per country
for (c in ddf$country){
  samples = rownames(metadata)[which(metadata$Country == c)]
  ddf[c,2] = length(samples)
}
ddf$samples_class = NA
ddf$samples_class[which(ddf$samples >=100)] = "100-1000"
ddf$samples_class[which(ddf$samples > 1000)] = ">1000"
ddf$samples_class[which(is.na(ddf$samples_class))] = "<100"
ddf$samples_class = factor(ddf$samples_class, levels=c("<100", "100-1000", ">1000"))

# edit world map template
data(wrld_simpl)
wrld_simpl@data$id <- wrld_simpl@data$NAME
wrld = fortify(wrld_simpl, region="id")
wrld = subset(wrld, id != "Antarctica")

# plot map
map.plot = ggplot() + 
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="white", color="grey40", linewidth=0.05) + 
  geom_map(data=ddf, map=wrld, aes(map_id=country, fill=samples_class),  color="grey40", linewidth=0.05, alpha=0.7) + 
  scale_fill_manual(values=c("darkolivegreen2", "darkgreen", "darkblue"), name="Number of samples") + 
  coord_map() + 
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        legend.position = "top")
ggsave("figures/map_metagenomes.pdf", height=5, width=8, dpi=300)