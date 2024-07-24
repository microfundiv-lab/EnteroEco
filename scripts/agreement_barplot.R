# load libraries
library(ggplot2)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/overlap/")
comb.coloniz.abund = read.csv("combined_coloniz-abund_all.csv", check.names=FALSE)
rownames(comb.coloniz.abund) = comb.coloniz.abund$Species_rep

comb.coloniz.abund.filt = read.csv("combined_coloniz-abund_filt.csv", check.names=FALSE)
rownames(comb.coloniz.abund.filt) = comb.coloniz.abund.filt$Species_rep

# check overlaps
tax = c("Entero", "Ecoli", "Kpneumo")
dataset = c("All", "Healthy-Adult")
disc.tax = c()
disc.dataset = c()

for (t in tax){
  for (d in dataset){
    tax.pos.dataset = rownames(comb.coloniz.abund)[which(rowSums(comb.coloniz.abund[,grep(t, colnames(comb.coloniz.abund))] > 0) > 0)]
    tax.neg.dataset = rownames(comb.coloniz.abund)[which(rowSums(comb.coloniz.abund[,grep(t, colnames(comb.coloniz.abund))] < 0) > 0)]
    inters = intersect(tax.pos.dataset, tax.neg.dataset)
    disc.dataset = c(disc.dataset, inters)
    
    tax.pos.taxon = rownames(comb.coloniz.abund)[which(rowSums(comb.coloniz.abund[,grep(d, colnames(comb.coloniz.abund))] > 0) > 0)]
    tax.neg.taxon = rownames(comb.coloniz.abund)[which(rowSums(comb.coloniz.abund[,grep(d, colnames(comb.coloniz.abund))] < 0) > 0)]
    inters = intersect(tax.pos.taxon, tax.neg.taxon)
    disc.tax = c(disc.tax, inters)
  }
}

disc.tax = unique(disc.tax)
disc.dataset = unique(disc.dataset)
disc.all = unique(c(disc.tax, disc.dataset))

# calculate disagreements per taxon
disagree.df = data.frame(matrix(NA, ncol=4, nrow=3))
colnames(disagree.df) = c("Agreement", "Disagree_Tax", "Disagree_Dataset", "Disagree_Both")
rownames(disagree.df) = c("Entero", "Ecoli", "Kpneumo")
for (n in 1:nrow(disagree.df)){
  tax = rownames(disagree.df)[n]
  tax.sig = rownames(comb.coloniz.abund)[which(rowSums(comb.coloniz.abund[,grep(tax, colnames(comb.coloniz.abund))] != 0) > 0)]
  tax.dis.dataset = intersect(tax.sig, disc.dataset)
  tax.dis.tax = intersect(tax.sig, disc.tax)
  tax.dis.both = intersect(tax.dis.dataset, tax.dis.tax)
  disagree.df$Disagree_Dataset[n] = length(which(!tax.dis.dataset %in% tax.dis.both))
  disagree.df$Disagree_Tax[n] = length(which(!tax.dis.tax %in% tax.dis.both))
  disagree.df$Disagree_Both[n] = length(tax.dis.both)
  disagree.df$Agreement[n] = length(tax.sig) - length(which(!tax.dis.dataset %in% tax.dis.both)) - length(which(!tax.dis.tax %in% tax.dis.both)) - length(tax.dis.both)
}

# plot results
disagree.melt = reshape2::melt(cbind(rownames(disagree.df), disagree.df))
colnames(disagree.melt)[1] = "Taxon"
bar.plot = ggplot(disagree.melt, aes(x=Taxon, y=value, fill=variable)) +
  geom_bar(stat="identity", alpha=0.6, width=0.6) +
  theme_classic() +
  scale_x_discrete(limits=c("Entero","Ecoli", "Kpneumo"), labels=c("Enterobacteriaceae", "E. coli", "K. pneumoniae")) +
  scale_fill_manual(values=c("darkgreen", "grey", "grey40", "black"), labels=c("Concordant", "Discordant (taxa)", "Discordant (dataset)", "Discordant (both)"), name="") +
  ylab("Number of species") +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.text = element_text(size=12), legend.position="right") +
  theme(axis.text.x = element_text(size=14, face="italic"))
ggsave("../figures/agreement_barplot.pdf", dpi=300, width=7, height=5)

# prop of agreement
prop.agree = (nrow(comb.coloniz.abund)-length(disc.all))/nrow(comb.coloniz.abund)*100
   
