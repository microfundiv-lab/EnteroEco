# load libraries
library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

# load metadata file
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
metadata = read.delim("metadata/metagenomes_07-2024_samples.tsv")
metadata = metadata[which(!is.na(metadata$Disease.name) & !is.na(metadata$Age.group) & !is.na(metadata$Continent)),]
rownames(metadata) = metadata$Sample

# metadata stats
state = reshape2::melt(table(metadata$Health.state))
state$Type = "Health state"

location = reshape2::melt(table(metadata$Continent))
location$Type = "Continent"

age = reshape2::melt(table(metadata$Age.group))
age$Type = "Age group"

stats = rbind(location, state, age)
stats = stats[which(!is.na(stats$Var1)),]
stats = stats[order(stats$Type,stats$value),]
stats$Var1 = factor(stats$Var1, levels=stats$Var1)
stats$value = stats$value/nrow(metadata)*100

# plot bargraph of metadata stats
fill_colors = c(brewer.pal(5,"Reds"), brewer.pal(3,"Greens")[2:3], brewer.pal(6,"Blues"))
age.plot = ggplot(stats[which(stats$Type == "Age group"),], aes(x=Type, y=value, fill=Var1)) +
  geom_bar(stat="identity", alpha=1) +
  coord_flip() +
  ylim(0,100.5) +
  theme_classic() +
  scale_fill_manual(values = brewer.pal(5,"Reds"), name="Age group") +
  ylab("% Metagenomic samples") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=14)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=12))

health.plot = ggplot(stats[which(stats$Type == "Health state"),], aes(x=Type, y=value, fill=Var1)) +
  geom_bar(stat="identity", alpha=1) +
  coord_flip() +
  ylim(0,100.5) +
  theme_classic() +
  scale_fill_manual(values = brewer.pal(3,"Greens")[2:3], name="Health state") +
  ylab("% Metagenomic samples") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=14)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=12))

cont.plot = ggplot(stats[which(stats$Type == "Continent"),], aes(x=Type, y=value, fill=Var1)) +
  geom_bar(stat="identity", alpha=1) +
  coord_flip() +
  ylim(0,100.5) +
  theme_classic() +
  scale_fill_manual(values = brewer.pal(6,"Blues"), name="Continent") +
  ylab("% Metagenomic samples") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=12))

# combine and save
ggarrange(age.plot, health.plot, cont.plot, nrow=3, align="v")
ggsave("figures/samples_distribution.pdf", height=5, width=7)






