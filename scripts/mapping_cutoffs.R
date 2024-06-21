# load libraries
library(data.table)
library(CoDaSeq)

# function
read_bwa = function (x) {
  bwa.data = as.data.frame(fread(x, check.names = FALSE))
  rownames(bwa.data) = bwa.data$Genome
  bwa.data = bwa.data[,metadata$Run]
  return(bwa.data)
}

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/bwa")
metadata = read.delim("../metadata/metagenomes_11-2023.tsv", check.names = TRUE)
rownames(metadata) = metadata$Run
bwa.cov.est = read_bwa("bwa_cov-est.csv")
bwa.cov.exp = read_bwa("bwa_cov-exp.csv")
bwa.counts = read_bwa("bwa_counts-filtered_batch-corr.csv")
outfile1 = "bwa_counts-filtered_batch-corr_test.csv"
outfile2 = "bwa_counts-filtered_batch-corr_clr.csv"

# define prevalence cut-offs
cov.est = 5 # Percentage
cov.ratio = 0.3 # 0-1
bwa.prev = bwa.cov.est
bwa.ratio = bwa.cov.est/bwa.cov.exp
bwa.prev[bwa.cov.est < cov.est | bwa.cov.est/bwa.cov.exp < cov.ratio]= 0
bwa.prev[bwa.prev != 0] = 1
bwa.counts = bwa.counts[rownames(bwa.prev), colnames(bwa.prev)]
bwa.counts[bwa.prev == 0] = 0

# clr transformation
bwa.clr = codaSeq.clr(bwa.counts + 0.5, samples.by.row=FALSE)

# save final tables
bwa.counts.fi = data.frame(cbind(rownames(bwa.counts), bwa.counts), check.names = FALSE)
colnames(bwa.counts.fi)[1] = "Genome"
write.table(bwa.counts.fi, file=outfile1, sep=",", quote=FALSE, row.names=FALSE)

bwa.clr.fi = data.frame(cbind(rownames(bwa.clr), bwa.clr), check.names = FALSE)
colnames(bwa.clr.fi)[1] = "Genome"
write.table(bwa.clr.fi, file=outfile2, sep=",", quote=FALSE, row.names=FALSE)