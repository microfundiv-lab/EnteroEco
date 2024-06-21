# load libraries
library(mikropml)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)

# function to get tpr and fpr per outcome
get_sens_spec_lookup = function(outcome, data){
  selected = data %>%
    select(!!as.name(outcome), observed)
  
  total = selected %>%
    count(observed) %>%
    pivot_wider(names_from=observed, values_from=n)
  
  selected %>%
    arrange(desc(!!as.name(outcome))) %>%
    mutate(is_outcome = (observed == outcome),
           tp = cumsum(is_outcome),
           fp = cumsum(!is_outcome),
           precision = tp / (tp + fp),
           sensitivity = tp / as.numeric(total[outcome]),
           fpr = fp / (sum(total)-as.numeric(total[outcome])),
           specificity = 1-fpr) %>% # TN/TN+FP = 1-FPR
    add_column(class = outcome) %>%
    select(precision, sensitivity, specificity, fpr, class)
}

# function to generate roc data per model
gen_roc_data = function(model, outcome, taxon, meta){
  select.model = readRDS(model)
  prob = predict(select.model$trained_model, select.model$test_data, type="prob")
  observed = select.model$test_data$Variable
  prob_obs = bind_cols(prob, observed=observed)
  roc.data = map_dfr(.x=outcome, .f=get_sens_spec_lookup, prob_obs)
  roc.data$Taxon = taxon
  roc.data$Dataset = meta
  return(roc.data)
}

# load best models
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/ml-species/best_models/")

# generate roc data
input.files = list.files(path = ".", pattern = ".Rds")
input.files = input.files[grep("HQ_",input.files)]
roc.list = lapply(input.files, function(x) {
  taxon = strsplit(x, "_")[[1]][2]
  dataset = strsplit(x, "_")[[1]][1]
  cat("Generating ROC data for", x, "...\n")
  roc.data = gen_roc_data(x, "Yes", taxon, dataset)
  return(roc.data)
})

roc.combined = as.data.frame(rbindlist(roc.list))
roc.combined$sensitivity = round(roc.combined$sensitivity, 2)
roc.combined$specificity = round(roc.combined$specificity, 2)
roc.agg.mean = aggregate(sensitivity ~ specificity + Taxon + Dataset, data=roc.combined, FUN=mean)
roc.agg.sd = aggregate(sensitivity ~ specificity + Taxon + Dataset, data=roc.combined, FUN=sd)
roc.agg.mean$sd = roc.agg.sd$sensitivity
roc.agg.mean$AUC = NA

# get AUCs
roc.auc.ecoli = read.csv("../HQ_Ecoli_results.csv")
roc.auc.ecoli = round(median(roc.auc.ecoli[which(roc.auc.ecoli$method == "xgbTree"),"AUC"]),3)
roc.agg.mean[which(roc.agg.mean$Taxon == "Ecoli"),"AUC"] = paste0("E. coli, AUC = ",roc.auc.ecoli)

roc.auc.entero = read.csv("../HQ_Entero_results.csv")
roc.auc.entero = round(median(roc.auc.entero[which(roc.auc.entero$method == "xgbTree"),"AUC"]),3)
roc.agg.mean[which(roc.agg.mean$Taxon == "Entero"),"AUC"] = paste0("Enterobacteriaceae, AUC = ",roc.auc.entero)

roc.auc.kp = read.csv("../HQ_Kpneumo_results.csv")
roc.auc.kp = round(median(roc.auc.kp[which(roc.auc.kp$method == "xgbTree"),"AUC"]),3)
roc.agg.mean[which(roc.agg.mean$Taxon == "Kpneumo"),"AUC"] = paste0("K. pneumoniae, AUC = ",roc.auc.kp)

# plot roc curve
roc.curve = ggplot(roc.agg.mean, aes(x=1-specificity, y=sensitivity, colour=AUC)) +
  geom_line(linewidth=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey") +
  theme_classic() + 
  theme(panel.grid.minor = element_blank()) + 
  ylim(0,1) +
  xlim(0,1) +
  ylab("True Positive Rate") + 
  xlab("False Positive Rate") +
  scale_colour_manual(values=c("steelblue", "darkgreen", "tomato")) +
  theme(legend.position="right") +
  guides(colour = guide_legend(override.aes = list(linewidth = 2))) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.background = element_blank(), strip.text = element_text(size=16)) +
  theme(legend.position=c(0.65,0.2), legend.box = "horizontal", legend.text=element_text(size=16),
        legend.title = element_blank()) +
  theme(axis.title.y = element_text(size=16)) + 
  theme(axis.text.y = element_text(size=14)) + 
  theme(axis.title.x = element_text(size=16)) + 
  theme(axis.text.x = element_text(size=14))
ggsave(filename="../../figures/ml-species_roc.pdf", width=7, height=5, dpi=300)
