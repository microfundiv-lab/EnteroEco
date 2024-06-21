# load libraries
library(reshape2)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# load data 
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/gem")
comb.abund = read.csv("../overlap/combined_coloniz-abund_filt.csv")
rownames(comb.abund) = comb.abund$Species_rep
comb.abund = comb.abund[which(comb.abund$Family != "Enterobacteriaceae"),]
comb.abund = comb.abund[,-1]

# add genome type
genome.metadata = read.delim("../metadata/genomes_uhgg-v1.2.tsv")[,c("Genome", "Genome_type")]
rownames(genome.metadata) = genome.metadata$Genome

# load gem data
gem.df = read.delim("cobrapy-M3constrain_results.tsv", header = F)
gem.df$V5 = gem.df$V2
gem.df$V2 = make.names(gem.df$V2)

# glm function
diff_metab = function(in_df, mtype) {
  gem.df = in_df[which(in_df$V4 == mtype),]
  gem.matrix = data.frame(acast(gem.df, V1 ~ V2, length))
  gem.matrix = gem.matrix[rownames(comb.abund),]
  gem.matrix[gem.matrix > 1] = 1
  metab.features = colnames(gem.matrix)
  metab2name = unique(gem.df[,c("V2", "V3", "V5")])
  rownames(metab2name) = metab2name$V2
  
  # prepare data
  final.df = merge(gem.matrix, comb.abund, by="row.names")
  rownames(final.df) = final.df$Row.names
  final.df = final.df[,-1]
  final.df$Genome_type = genome.metadata[rownames(final.df),"Genome_type"]
  
  # run glm
  glm.summary = lapply(metab.features, function(x) {
    cat(paste("Running glm for", x, "\n", sep=" "))
    glm.formula = formula(paste(paste0(x, " ~ Classification + Genome_type")))
    glm.out = summary(glm(glm.formula, data = final.df, family = "binomial"))
    return(glm.out$coefficients["ClassificationCo-excluder",-3])
  })
  
  # process output
  raw.output = data.frame(t(data.frame(glm.summary)))
  rownames(raw.output) = metab.features
  raw.output$FDR = p.adjust(raw.output[,3])
  sign.output = raw.output[which(raw.output$FDR < 0.05),]
  sign.output = sign.output[order(sign.output$FDR),]
  sign.fi = cbind(rownames(sign.output), sign.output)
  sign.fi$Estimate = -sign.fi$Estimate
  colnames(sign.fi) = c("Metabolite", "Estimate", "Std.Error", "Pvalue", "FDR")
  sign.fi = sign.fi[order(sign.fi$Estimate, decreasing=TRUE),]
  sign.fi$Original = metab2name[sign.fi$Metabolite, "V5"]
  return(sign.fi)
}

# heatmap function
plot_heatmap = function(in_df, mtype, glm_df) {
  gem.df = in_df[which(in_df$V4 == mtype),]
  gem.matrix = data.frame(acast(gem.df, V1 ~ V2, length))
  gem.matrix = gem.matrix[rownames(comb.abund),]
  gem.matrix[gem.matrix > 1] = 1
  metab.features = colnames(gem.matrix)
  metab2name = unique(gem.df[,c("V2", "V3")])
  rownames(metab2name) = metab2name$V2
  
  # prepare data
  final.df = merge(gem.matrix, comb.abund, by="row.names")
  rownames(final.df) = final.df$Row.names
  final.df = final.df[,-1]
  final.df$Genome_type = genome.metadata[rownames(final.df),"Genome_type"]
  
  # prepare colors
  neg = rownames(final.df)[which(final.df$Classification == "Co-excluder")]
  pos = rownames(final.df)[which(final.df$Classification == "Co-colonizer")]
  annot = final.df[c(neg, pos), c("Classification", "Genome_type", "Order"), drop=FALSE]
  annot$Classification = gsub("Positive", "Co-colonizer", annot$Classification)
  annot$Classification = gsub("Negative", "Co-excluder", annot$Classification)
  class.colors = c("steelblue", "tomato")
  names(class.colors) = c("Co-excluder", "Co-colonizer")
  status.colors = c("steelblue1", "palegreen2")
  names(status.colors) = c("Isolate", "MAG")
  top.taxa = names(sort(table(annot$Order), decreasing = TRUE)[1:5])
  annot$Order = ifelse(annot$Order %in% top.taxa, annot$Order, "Other")
  top.taxa.colors = c(brewer.pal(length(top.taxa), "Set3"), "lightgrey")
  names(top.taxa.colors) = c(top.taxa, "Other")
  annot.colors = list(Classification=class.colors, Genome_type=status.colors, Order=top.taxa.colors)
  gem.annot = data.frame(ifelse(glm_df$Estimate > 0, "Co-colonizer", "Co-excluder"))
  rownames(gem.annot) = rownames(glm_df)
  colnames(gem.annot) = "Classification"
  
  # heatmap plot
  sign.df = t(final.df[c(neg, pos), glm_df$Metabolite])
  pheatmap(sign.df, show_colnames = FALSE, show_rownames = TRUE, cluster_rows = TRUE,
           annotation_col = annot, annotation_colors = annot.colors,
           annotation_names_row = FALSE, annotation_names_col=TRUE,
           color=c("lightgrey", "darkblue"), labels_row = metab2name[rownames(sign.df),"V3"],
           filename=paste0("../figures/gem-heatmap_", mtype, ".pdf"), width = 15, height=13)
}

# run functions
uptake.df = diff_metab(gem.df, "Uptake")
secretion.df = diff_metab(gem.df, "Secretion")
plot_heatmap(gem.df, "Uptake", uptake.df)
plot_heatmap(gem.df, "Secretion", secretion.df)

# check overlap
overlap.df = merge(uptake.df, secretion.df, by="Metabolite")

# save list of metabolites
all.metab = unique(c(uptake.df$Original, secretion.df$Original))
write.table(all.metab, file="metab_candidates.txt", col.names=FALSE, quote=FALSE, row.names=FALSE)
