# load libraries
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(tidyr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
all.results.raw = read.delim("timeseries/maaslin2/all_results.tsv")
all.results.raw$coef = -all.results.raw$coef

# get taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
species.data = read.delim("metadata/species_uhgg_v1.2.tsv")
species.data = unique(species.data[,c("Species_rep", "Lineage")])
species.data = separate(data = species.data, col = Lineage, sep = ";", into = ranks)
rownames(species.data) = species.data$Species_rep
species.data = data.frame(lapply(species.data, function(x) { gsub(".__", "", x)}))
rownames(species.data) = species.data$Species_rep

# color dictionary
families = sort(c("Lachnospiraceae", "Ruminococcaceae", "Acutalibacteraceae", "Enterococcaceae", "Enterobacteriaceae", 
             "Bacteroidaceae", "UBA1381", "Butyricicoccaceae","Veillonellaceae"), decreasing=TRUE)
colors = brewer.pal(length(families), "Set3")
names(colors) = families

# function to process results
process_maaslin2 = function(status, label) {
  all.results = all.results.raw[which(all.results.raw$value == status),]
  all.results$FDR = p.adjust(all.results$pval, method="fdr")
  ori.species = scan("metadata/species_list.txt", what="")
  all.results = all.results[which(all.results$feature %in% ori.species),]
  sig.results = all.results[which(all.results$FDR < 0.05),]
  
  # compare results
  timeseries.cocol = sig.results[which(sig.results$coef > 0), "feature"]
  timeseries.coexcl = sig.results[which(sig.results$coef < 0), "feature"]
  ori.results = read.csv("overlap/combined_coloniz-abund_all.csv")
  rownames(ori.results) = ori.results$Species_rep
  ori.cocol = ori.results[which(ori.results$Classification == "Co-colonizer"), "Species_rep"]
  ori.coexcl = ori.results[which(ori.results$Classification == "Co-excluder"), "Species_rep"]
  
  # report overlaps
  agree.coexcl = unique(intersect(timeseries.coexcl, ori.coexcl))
  agree.cocol = unique(intersect(timeseries.cocol, ori.cocol))
  disagree.coexcl = unique(intersect(timeseries.cocol, ori.coexcl))
  disagree.cocol = unique(intersect(timeseries.coexcl, ori.cocol))
  missing = length(unique(sig.results[which(!sig.results$feature %in% ori.results$Species_rep),"feature"]))

  overlap.df = data.frame(matrix(nrow=3, ncol=3))
  colnames(overlap.df) = c("Ref", "Species", "Number")
  overlap.df$Ref = "Significant species"
  overlap.df$Species = c("Discordant", "Concordant", "Missing")
  overlap.df$Species = factor(overlap.df$Species, levels=c("Discordant", "Missing", "Concordant"))
  overlap.df$Number[which(overlap.df$Species == "Discordant")] = length(c(disagree.coexcl, disagree.cocol))
  overlap.df$Number[which(overlap.df$Species == "Concordant")] = length(c(agree.coexcl, agree.cocol))
  overlap.df$Number[which(overlap.df$Species == "Missing")] = missing
  overlap.df$Percentage = paste0(round(overlap.df$Number/sum(overlap.df$Number)*100, digits=0),"%")
  selected.species = unique(c(agree.coexcl, agree.cocol))
  
  # chisq test
  sig.ori = length(c(agree.coexcl, agree.cocol))
  sig.times = length(unique(sig.results$feature))-sig.ori
  eligible.species = intersect(ori.results$Species_rep, unique(all.results$feature))
  nonsig.ori = length(eligible.species)-sig.ori
  nonsig.times = length(unique(all.results$feature))-nonsig.ori-sig.times-sig.ori
  chisq = chisq.test(rbind(c(sig.ori, sig.times), c(nonsig.ori, nonsig.times)))
  print(chisq)
  
  barplot.df = all.results.raw[which(all.results.raw$feature %in% selected.species & all.results.raw$value == status),]
  rownames(barplot.df) = barplot.df$feature
  barplot.df$Family = species.data[rownames(barplot.df), "Family"]
  barplot.df$Species = species.data[rownames(barplot.df), "Species"]
  barplot.df$Species = gsub("_", " ", barplot.df$Species)
  barplot.df = barplot.df[order(barplot.df$coef, decreasing=FALSE),]
  
  # plots
  species.plots = ggplot(barplot.df, aes(x=Species, y=log2(exp(coef)), fill=Family)) +
    geom_bar(stat="identity", alpha=0.7) +
    geom_errorbar(aes(ymin = log2(exp(coef))-log2(exp(stderr)), ymax=log2(exp(coef))+log2(exp(stderr))), width=0.1) +
    coord_flip() +
    theme_classic() +
    ylab(bquote("Effect size (log"[2]*" fold-change)")) +
    scale_fill_manual(values = colors) +
    scale_x_discrete(limits = barplot.df$Species) +
    theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size = 12, face = 'italic'), axis.title.y = element_blank(), axis.title.x = element_text(size=14),
          legend.text = element_text(size = 12, face = 'italic', hjust = 0), legend.title = element_text(size=12), legend.position = "right")
  
  overlap.plot = ggplot(overlap.df, aes(x=Ref, y=Number, fill=Species, label=Percentage)) +
    geom_bar(stat="identity", alpha=0.5) +
    coord_flip() +
    theme_classic() +
    ggtitle(label) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values=c("tomato", "grey", "darkgreen"), labels=c("Discordant", "Missing", "Concordant"), name="") +
    ylab("Number of species") +
    theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size = 12), axis.title.y = element_blank(), axis.title.x = element_text(size=14),
          legend.text = element_text(size = 12), legend.title = element_text(size=12), legend.position = "right")
  ggarrange(overlap.plot, species.plots, nrow=2, heights=c(0.7,2))
}

# generate tables
overlap.control = process_maaslin2("CPE_negative_control", "CPE-negative (controls)")
overlap.index = process_maaslin2("CPE_negative_index", "CPE-negative (index)")
ggarrange(overlap.control, overlap.index, ncol=2)
ggsave(filename = "figures/timeseries_validation.pdf", dpi=300, width=14, height=6)
