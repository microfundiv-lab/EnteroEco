# load libraries
library(ggplot2)
library(ape)
library(igraph)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data")
metadata = read.delim("metadata/metagenomes_11-2023.tsv")
rownames(metadata) = metadata$Run

# load mlst results
st_results = read.delim("metamlst/escherichia1_report.txt")

# confidence filter
st_results = st_results[which(st_results$Confidence > 90),]

# parse results
st_results$Sample = metadata[st_results$X,"Sample"]
st_results = unique(st_results[,c("ST", "Sample")])
st_results$Class = ifelse(as.numeric(st_results$ST) > 100000, "Novel", "Known")
st_alleles = read.delim("metamlst/escherichia1_ST.txt")
st_alleles = st_alleles[which(st_alleles$ST %in% st_results$ST),]
rownames(st_alleles) = st_alleles$ST
st_alleles = st_alleles[,-1]

# st frequencies
st_freq = data.frame(table(st_results$ST), stringsAsFactors = FALSE)
st_freq$Var1 = as.character(st_freq$Var1)
rownames(st_freq) = st_freq$Var1
st_freq = st_freq[rownames(st_alleles),]
st_freq$Class = ifelse(as.numeric(st_freq$Var1) > 100000, "Novel", "Known")

# calculate distances
d.mlst.distances <- matrix(0, ncol=nrow(st_alleles), nrow=nrow(st_alleles))
for (i in 1:(nrow(st_alleles)-1)) {
  for (j in (i+1):nrow(st_alleles)){
    d.mlst.distances[i,j] = sum(st_alleles[i, ] != st_alleles[j, ])
    d.mlst.distances[j,i] = sum(st_alleles[i, ] != st_alleles[j, ])
  }
}
d.mlst.distances.euclidean = dist(d.mlst.distances)

# get minimum spanning tree and convert to igraph graph object
d.mst = ape::mst(d.mlst.distances.euclidean)
g.mst = graph_from_adjacency_matrix(d.mst, mode=c('undirected'))

# Set vertex size to ST frequency
V(g.mst)$size = log(round( 20 * (st_freq$Freq / max(st_freq$Freq)) ) + 1)*2
V(g.mst)$label = ifelse(st_freq$Freq > 5, st_freq$Var1, "")
V(g.mst)$color = ifelse(st_freq$Class == "Known", "steelblue1", "darkgreen")

# plot graph
set.seed(3000)
pdf("figures/mlst_graph.pdf", height=15, width=7)
plot.igraph(g.mst, edge.arrow.size=0.5, vertex.frame.color="gray", 
            vertex.label.family = "sans", vertex.label.color = "black", vertex.label.dist = 0.75, vertex.label.degree = -pi/2,
            layout=layout.fruchterman.reingold, margin=0)
dev.off()

# check stats
st_top10 = st_freq[order(st_freq$Freq, decreasing=TRUE),][1:10,]
props = table(st_results$Class)/nrow(st_results)*100
