#!/usr/bin/env Rscript

# load libraries
library(NetCoMi)

# load data
cat("Loading data ...\n")
setwd("/rds/project/rds-aFEMMKDjWlo/aalmeida/project_ecoevo/netcomi")
entero_no = read.delim("healthy-adult_entero-neg_corr-filtered.tsv", header = T)
rownames(entero_no) = entero_no$X.OTU.ID
entero_no = entero_no[,-1]
entero_yes = read.delim("healthy-adult_entero-pres_corr-filtered.tsv", header = T)
rownames(entero_yes) = entero_yes$X.OTU.ID
entero_yes = entero_yes[,-1]

# filter to common taxa
common.taxa = intersect(colnames(entero_no), colnames(entero_yes))
entero_no = entero_no[common.taxa, common.taxa]
entero_yes = entero_yes[common.taxa, common.taxa]

# construct networks
cat("Running netConstruct ...\n")
net.construct = netConstruct(data = as.matrix(entero_no), data2 = as.matrix(entero_yes),
                          dataType = "correlation", sparsMethod = "none", seed = 1234)

# identify differentially associated taxa
cat("Performing differential analysis ...\n")
diff.res = diffnet(x=net.construct, diffMethod="fisherTest", adjust = "lfdr", n1=1960, n2=3145)
diffmat_sums = rowSums(diff.res$diffAdjustMat)
diff_asso_names = names(diffmat_sums[diffmat_sums > 0])

# analyse networks
cat("Analysing networks ...\n")
net.analyse = netAnalyze(net.construct, centrLCC = FALSE, avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE, clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "eigenvector"), hubQuant = 0.9,
                           lnormFit = TRUE, normDeg = FALSE, normBetw = FALSE,
                           normClose = FALSE, normEigen = FALSE)

# compare networks
cat("Comparing networks ...\n")
net.compare = netCompare(net.analyse)
summary(net.compare)

# visualize networks
cat("Building network plot ...\n")
png("netcomi_networks.png", res=300, width = 12, height = 6, units="in")
plot(net.analyse, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", groupNames = c("",""), negDiffCol = FALSE, 
     edgeWidth=10, cut=0, edgeTranspLow=0, edgeFilter = "threshold", edgeFilterPar = 0.1, nodeColor = "darkgreen",
     nodeFilter = "names", nodeFilterPar = diff_asso_names, nodeSize = "degree", posCol = "darkgrey", labels = FALSE, highlightHubs = FALSE, cexNodes = 1.5, hubBorderCol  = "gray40")
dev.off()

# save environment
cat("Saving final environment ...\n")
save.image(file='netcomi_entero.RData')
