# load libraries
library(ape)
library(vegan)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
adonis_test = function(strain) {
  tree.file = read.tree(paste0("itol/strains/", strain, ".nwk"))
  cophenetic.data = cophenetic.phylo(tree.file)
  gene.counts = read.table(paste0("strains/", strain, "_uhgp-filt.tsv"))
  gene.counts.df = as.data.frame(table(gene.counts$V2))
  rownames(gene.counts.df) = gene.counts.df$Var1
  cophenetic.subs = cophenetic.data[rownames(gene.counts.df),rownames(gene.counts.df)]
  
  # perform adonis
  adonis.res = adonis2(cophenetic.subs ~ Freq, data=gene.counts.df)
  return(adonis.res)
}
gut718 = adonis_test("GUT_GENOME000718")
gut962 = adonis_test("GUT_GENOME096210")
