# load libraries
library(reshape2)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# load data 
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/genofan/")
comb.abund = read.csv("../overlap/combined_coloniz-abund_filt.csv")
rownames(comb.abund) = comb.abund$Species_rep
comb.abund = comb.abund[which(comb.abund$Family != "Enterobacteriaceae"),]
comb.abund = comb.abund[,-1]

# add genome type
genome.metadata = read.delim("../metadata/genomes_uhgg-v1.2.tsv")[,c("Genome", "Genome_type")]
rownames(genome.metadata) = genome.metadata$Genome

# parse kegg results
ko.classes.ori = read.delim("KO_Orthology_ko00001.txt", header=FALSE)
ko.classes = separate(ko.classes.ori, V4, "V5", sep=" ")
ko.desc = cbind(ko.classes.ori, ko.classes)
ko.desc = unique(ko.desc[,c("V5", "V4")])
rownames(ko.desc) = ko.desc$V5

# load kegg data
keggor.df = read.delim("kegg_orthologs.tsv", sep="", header = F)
keggor.df$V1 = gsub("_\\d+$", "", keggor.df$V1)
keggor.matrix = data.frame(acast(keggor.df, V2 ~ V1, length))
keggor.matrix = keggor.matrix[,rownames(comb.abund)]
keggor.matrix = t(keggor.matrix[which(rowSums(keggor.matrix > 0) > ncol(keggor.matrix)*0.01),])
keggor.matrix[keggor.matrix > 0] = 1
kegg.features = colnames(keggor.matrix)

# prepare data
final.df = merge(keggor.matrix, comb.abund, by="row.names")
rownames(final.df) = final.df$Row.names
final.df = final.df[,-1]
final.df$Genome_type = genome.metadata[rownames(final.df),"Genome_type"]

# run glm
glm.summary = lapply(kegg.features, function(x) {
  cat(paste("Running glm for", x, "\n", sep=" "))
  glm.formula = formula(paste(paste0(x, " ~ Classification + Genome_type")))
  glm.out = summary(glm(glm.formula, data = final.df, family = "binomial"))
  return(glm.out$coefficients["ClassificationCo-excluder",-3])
})

# process output
raw.output = data.frame(t(data.frame(glm.summary)))
rownames(raw.output) = kegg.features
raw.output$FDR = p.adjust(raw.output[,3])
sign.output = raw.output[which(raw.output$FDR < 0.05),]
sign.output = sign.output[order(sign.output$FDR),]
sign.fi = cbind(rownames(sign.output), sign.output)
sign.fi$Estimate = -sign.fi$Estimate
colnames(sign.fi) = c("KEGG_Ortholog", "Estimate", "Std.Error", "Pvalue", "FDR")
sign.fi$Description = ko.desc[sign.fi$KEGG_Ortholog,"V4"]
write.table(sign.fi, file="ko-glm_results.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# select top ko and parse description
ko.pos = rownames(sign.fi)[which(sign.fi$Estimate > 0)]
ko.pos.top = sign.fi[order(sign.fi$Estimate, decreasing=TRUE)[1:20],]
ko.neg = rownames(sign.fi)[which(sign.fi$Estimate < 0)]
ko.neg.top = sign.fi[order(sign.fi$Estimate, decreasing=FALSE)[1:20],]
ko.top = rbind(ko.pos.top, ko.neg.top)
ko.top = ko.top[order(ko.top$Estimate, decreasing=TRUE),]
ko.top$Description = ko.desc[ko.top$KEGG_Ortholog,"V4"]
ko.top$Description = gsub(".*;","",ko.top$Description)
ko.top$Description = gsub("\\[.*","",ko.top$Description)
ko.top$Description = str_to_title(ko.top$Description)
ko.top$Description = gsub("16s","16S",ko.top$Description)
ko.top$Description = gsub("30s","30S",ko.top$Description)
ko.top$Description = gsub("Trna","tRNA",ko.top$Description)
ko.top$Description = gsub("Rna","RNA",ko.top$Description)
ko.top$Description = gsub("Mfs","MFS",ko.top$Description)
ko.top$Description = gsub("Udp","UDP",ko.top$Description)

# heatmap colors
neg = rownames(final.df)[which(final.df$Classification == "Co-excluder")]
pos = rownames(final.df)[which(final.df$Classification == "Co-colonizer")]
annot = final.df[c(neg, pos), c("Classification", "Genome_type", "Order"), drop=FALSE]
class.colors = c("steelblue", "tomato")
names(class.colors) = c("Co-excluder", "Co-colonizer")
status.colors = c("steelblue1", "palegreen2")
names(status.colors) = c("Isolate", "MAG")
top.taxa = names(sort(table(annot$Order), decreasing = TRUE)[1:5])
annot$Order = ifelse(annot$Order %in% top.taxa, annot$Order, "Other")
top.taxa.colors = c(brewer.pal(length(top.taxa), "Set3"), "lightgrey")
names(top.taxa.colors) = c(top.taxa, "Other")
annot.colors = list(Classification=class.colors, Genome_type=status.colors, Order=top.taxa.colors)
kegg.annot = data.frame(ifelse(ko.top$Estimate > 0, "Co-colonizer", "Co-excluder"))
rownames(kegg.annot) = rownames(ko.top)
colnames(kegg.annot) = "Classification"

# heatmap plot
sign.df = t(final.df[c(neg, pos), ko.top$KEGG_Ortholog])
pheatmap(sign.df, show_colnames = FALSE, show_rownames = TRUE, cluster_rows = TRUE,
         annotation_col = annot, annotation_colors = annot.colors,
         annotation_names_row = FALSE, annotation_names_col=TRUE,
         color=c("lightgrey", "darkgreen"), labels_row = ko.top$Description,
         filename="../figures/ko-heatmap_top20.pdf", width = 15, height=7)

# investigate categories
ko.pos.categories = ko.classes[which(ko.classes$V5 %in% ko.pos.top$KEGG_Ortholog),]
ko.neg.categories = ko.classes[which(ko.classes$V5 %in% ko.neg.top$KEGG_Ortholog),]
