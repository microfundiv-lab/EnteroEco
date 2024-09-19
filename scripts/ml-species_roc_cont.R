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
gen_roc_data = function(model, outcome, continent){
  select.model = readRDS(model)
  prob = predict(select.model$trained_model, select.model$test_data, type="prob")
  observed = select.model$test_data$Variable
  prob_obs = bind_cols(prob, observed=observed)
  roc.data = map_dfr(.x=outcome, .f=get_sens_spec_lookup, prob_obs)
  roc.data$Continent = continent
  return(roc.data)
}

# load best models
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_QiYin_Antipath/data/ml-species/continent_within/")

# generate roc data
input.files = list.files(path = ".", pattern = ".Rds", recursive = TRUE)
roc.list = lapply(input.files, function(x) {
  cat("Generating ROC data for", x, "...\n")
  continent = gsub("_.*", "", x)
  roc.data = gen_roc_data(x, "Yes", continent)
  roc.data$TPR = 1-roc.data$specificity
  return(roc.data)
})

roc.combined = as.data.frame(rbindlist(roc.list))
roc.combined$TPR = round(roc.combined$TPR, digits=3)
roc.combined$sensitivity = round(roc.combined$sensitivity, digits=3)
roc.agg.mean = aggregate(sensitivity ~ TPR + Continent, data=roc.combined, FUN=mean)
roc.agg.mean$AUC = NA

# ensure monotonicity
for (cont in unique(roc.agg.mean$Continent)) {
  cont.rows = which(roc.agg.mean$Continent == cont)
  for (n in cont.rows[2:length(cont.rows)]){
    roc.agg.mean[n,"sensitivity"] = ifelse(roc.agg.mean[n-1,"sensitivity"] > roc.agg.mean[n,"sensitivity"], roc.agg.mean[n-1,"sensitivity"], roc.agg.mean[n,"sensitivity"])
  }
}

# get AUCs
perf.files = list.files(path = ".", pattern = "performance_results.csv", recursive=TRUE)
for (perf in perf.files) {
  roc.auc.cont = read.csv(perf)
  continent = gsub("_.*", "", perf)
  cont.space = gsub("NorthAmerica", "North America", continent)
  cont.space = gsub("SouthAmerica", "South America", cont.space)
  auc.value = round(median(roc.auc.cont[,"AUC"], na.rm=TRUE),3)
  roc.agg.mean[which(roc.agg.mean$Continent == continent),"AUC"] = paste0(cont.space, ", AUROC = ", auc.value)
}

# plot roc curve
roc.curve = ggplot(roc.agg.mean, aes(x=TPR, y=sensitivity, colour=AUC)) +
  geom_line(linewidth=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey") +
  theme_classic() + 
  theme(panel.grid.minor = element_blank()) + 
  ylim(0,1) +
  xlim(0,1) +
  ylab("True Positive Rate") + 
  xlab("False Positive Rate") +
  scale_colour_manual(values=brewer.pal(10, "Set3")) +
  theme(legend.position="right") +
  guides(colour = guide_legend(override.aes = list(linewidth = 2))) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.background = element_blank(), strip.text = element_text(size=16)) +
  theme(legend.position=c(0.65,0.2), legend.box = "horizontal", legend.text=element_text(size=14),
        legend.title = element_blank()) +
  theme(axis.title.y = element_text(size=16)) + 
  theme(axis.text.y = element_text(size=14)) + 
  theme(axis.title.x = element_text(size=16)) + 
  theme(axis.text.x = element_text(size=14))
ggsave(filename="../../figures/ml-species_roc_cont.pdf", width=7, height=5, dpi=300)
