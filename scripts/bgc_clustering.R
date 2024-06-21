# load libraries
library(dplyr)
library(readr)
library(data.table)
library(igraph)
library(Matrix)
library(reshape2)
library(ggplot2)
library(GGally)
library(RColorBrewer)
library(stringr)
library(ggpubr)
library(matrixStats)
library(ggrastr)
library(ggtext)

# load metadata and bgc blast files
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/")
metadata = read.csv("overlap/combined_coloniz-abund_filt.csv")
blastn.all = read.delim("genofan/all_vs_all_ani.tsv", header=TRUE)
blastn.all = blastn.all[which(blastn.all$pid >= 50),]
blastn.all$cov = rowMaxs(as.matrix(blastn.all[,c("qcov", "tcov")]))

# select BGC (cyclic lactone)
qname.keep = blastn.all$qname[which(grepl("cyclic-lactone-autoinducer", blastn.all$qname))]
tname.keep = blastn.all$tname[which(grepl("cyclic-lactone-autoinducer", blastn.all$tname))]
blastn.all = blastn.all[which(blastn.all$qname %in% qname.keep | blastn.all$tname %in% tname.keep),]

# check cov distribution
hist.plot = ggplot(blastn.all, aes(x=cov)) +
  geom_density(fill="darkgrey") +
  theme_classic() +
  xlab("Coverage (%)") +
  ylab("Density") +
  geom_vline(xintercept=50, linetype="dashed", linewidth=0.1) +
  theme(axis.title = element_text(size=14)) +
  theme(axis.text = element_text(size=12))
ggsave(filename = "figures/bgc_cluster-thresh.pdf", height=6, width=7, dpi=300)

# make matrix symmetrical
blastn.symm = unique(rbind(blastn.all[,c("qname", "tname", "cov")], rbind(blastn.all[,c("tname", "qname", "cov")])))

# convert to distance and filter
blastn.matrix = as.data.frame(acast(qname ~ tname, by="cov", data=blastn.symm))
blastn.matrix[is.na(blastn.matrix)] = 0
blastn.matrix[blastn.matrix < 50] = 0

# convert to graph
order.bgcs = colnames(blastn.matrix)
blastn.matrix = blastn.matrix[order.bgcs, order.bgcs]
net = graph.adjacency(as.matrix(blastn.matrix), weighted=TRUE, diag = FALSE)

# remove vertices absent from candidate list
all.genomes = gsub("(GUT_GENOME\\d+)_.*", "\\1", V(net)$name)
to.remove = which(!all.genomes %in% metadata$Species_rep | !grepl("cyclic-lactone-autoinducer", V(net)$name))
net.filt = delete.vertices(net, to.remove)

# define graph metadata
genome.names = gsub("(GUT_GENOME\\d+)_.*", "\\1", V(net.filt)$name)
V(net.filt)$Family = metadata[match(genome.names, metadata$Species_rep),"Family"]
V(net.filt)$Classification = metadata[match(genome.names, metadata$Species_rep),"Classification"]
V(net.filt)$Classification = gsub("Negative", "Co-excluder", V(net.filt)$Classification)
V(net.filt)$Classification = gsub("Positive", "Co-colonizer", V(net.filt)$Classification)
class.color = c("Co-excluder" = "steelblue", "Co-colonizer" = "tomato")
class.shape = c("Co-excluder" = 21, "Co-colonizer" = 25)
family.color = rev(brewer.pal(length(unique(V(net.filt)$Family)), "Set3"))
names(family.color) = unique(V(net.filt)$Family)
family.color = sort(family.color)

# plot graph
set.seed(100)
net.plot = ggnet2(net.filt, edge.color = "darkgrey", edge.alpha = 0.8, edge.size=abs(E(net.filt)$weight)/50) +
  geom_point(aes(fill=V(net.filt)$Family, shape=V(net.filt)$Classification), size=5, colour="lightgrey") +
  scale_fill_manual(values=family.color, name="Family") +
  scale_shape_manual(values=class.shape, name="Classification") +
  guides(fill = guide_legend(override.aes = list(size=5, shape=21)),
        shape = guide_legend(override.aes = list(size=5))) +
  theme(legend.text = element_text(size=12), legend.title = element_text(size=12))
ggsave(filename = "figures/bgc_network.pdf", height=6, width=8, bg="white", dpi=300)

# load mibig results
mibig.df = read.delim("genofan/mibig_result.tsv", header = F)
mibig.df = mibig.df[,c("V1", "V2", "V3")]
colnames(mibig.df) = c("query", "known.bgc", "identity")
mibig.df$locus = gsub(".*?(GUT_GENOME[0-9]+_[0-9]+).*", "\\1", mibig.df$query)
locus.df = read.delim("genofan/bgc_region2locus.tsv", header = F)
colnames(locus.df) = c("region", "locus")
target.df = merge(mibig.df, locus.df, by="locus")

# get cluster memberships
node.cluster = as.data.frame(components(net.filt)$membership)
colnames(node.cluster) = "cluster.group"
node.cluster$region = gsub("^(.*?)__.*", "\\1", names(components(net.filt)$membership))
target.df = merge(target.df, node.cluster, by = "region")
target.df$Species_rep = gsub("(GUT_GENOME\\d+)_.*", "\\1", target.df$locus)
target.df$Family = metadata[match(target.df$Species_rep , metadata$Species_rep),"Family"]

# plot distribution
dis.mibig = ggplot(target.df, aes(x = reorder(cluster.group, identity, FUN = median), y = identity, fill = Family)) +
  geom_point_rast(alpha=0.1, size=0.2, position = position_jitter(width = 0.2)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  theme_classic() +
  coord_flip() +
  ylab("Amino acid identity (%)") +
  xlab("BGC family") +
  scale_fill_manual(values=family.color) +
  guides(fill="none") +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=10)) + 
  theme(axis.title.x = element_text(size=16)) +
  theme(axis.text.x = element_text(size=14))
ggsave(filename = "figures/bgc_boxplot.pdf", height=10, width=4, bg="white", dpi=300)

# get numbers
auto.bgcs.all = length(V(net.filt))
auto.bgcs.clsts = components(net.filt)$no
top.families = sort(table(V(net.filt)$Family))
top.families.prop = sort(table(V(net.filt)$Family)/auto.bgcs.all*100)
ani.iqr = quantile(target.df$identity)