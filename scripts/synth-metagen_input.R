# load libraries
library(data.table)

# do not use scientific notation
options(scipen=999)

# load input
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
abund.data = as.data.frame(fread("bwa/bwa_counts-filtered_samples.csv"))
rownames(abund.data) = abund.data$Genome
abund.data = abund.data[,-1]
genome.sizes.all = read.delim("metadata/species_uhgg_v1.2.tsv")
rownames(genome.sizes.all) = genome.sizes.all$Species_rep
genome.sizes.noentero = genome.sizes.all[which(!grepl("f__Enterobacteriaceae", genome.sizes.all$Lineage)),]

# select top taxa
n = 50
top.taxa = names(sort(rowSums(abund.data[rownames(genome.sizes.noentero),] > 0), decreasing=TRUE)[1:n])
synth.abund = genome.sizes.noentero[top.taxa,c("Species_rep", "Length")]

# function to add one Enterobacteriaceae species
input.depths = c(12639672, 30869056, 50929724) # read depth quartiles (25%, 50%, 75%)
gen_synth_comm = function(input.species, input.abund) {
  synth.abund.new = rbind(c(input.species, genome.sizes.all[input.species,"Length"]),synth.abund)
  rownames(synth.abund.new)[1] = input.species
  synth.abund.new$Abundance = (1-(input.abund/100))/n
  synth.abund.new$Abundance[1] = input.abund/100
  synth.abund.new$Read.counts = round(as.numeric(synth.abund.new$Length)*synth.abund.new$Abundance, digits=0)
  for (depth in input.depths) {
    depth.ratio = (depth/2)/sum(synth.abund.new$Read.counts)
    synth.abund.new[,paste(input.species, input.abund, paste(round(depth/1000, digits=0), "k", sep=""), sep="_")] = round(synth.abund.new$Read.counts*depth.ratio, digits=0)
  }
  synth.melt = reshape2::melt(synth.abund.new[,which(!colnames(synth.abund.new) %in% c("Length", "Abundance", "Read.counts"))])
  write.table(synth.melt, file=paste0("synth_metagen/", paste(input.species, input.abund, sep="_"), ".tsv"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}

# run for selected species and abundances
entero.selected = c("GUT_GENOME144544", "GUT_GENOME147598","GUT_GENOME143760", "GUT_GENOME000563", "GUT_GENOME231247")
abund.selected = c(1, 0.5, 0.1, 0.05, 0.01, 0.0075, 0.005, 0.004, 0.003, 0.002, 0.001)
for (species in entero.selected) {
  for (abund in abund.selected) {
    gen_synth_comm(species, abund)
  }
}