library("ggsci")
setwd('/Users/Yuan/Desktop/Variant_calling/')

library(tidyverse)
library(ggplot2)
library(viridis)
library(reshape2)


### imputed vs. original called
readin_imputed_original_inconsistent <- function(pattern){
  files = list.files('data/performance/')
  files = files[grepl("variants_numbers", files)]
  files = files[grepl(pattern, files)]
  dat = NULL
  for(f in files){
    dati = read.table(paste0('data/performance/', f), sep='\t', header = T, stringsAsFactors = F)
    sample = strsplit(f, '_')[[1]][1]
    dati$sample = sample
    dat = rbind(dat, dati)
  }
  dat$group = pattern
  return(dat)
}


### Number of recovered variants
pattern = 'minDP4'
dat = readin_imputed_original_inconsistent(pattern)
for(col in colnames(dat)[1:5]){
  dat[,col] = dat[,col] / dat$Golden_standard
}
number_recovered = dat[, c("sample", "group", c("Only_from_gc", "Two_consistent", "Two_inconsistent", "Only_from_imputation"))]
number_recovered = melt(number_recovered, id.vars = c("sample", "group"))
number_recovered$variable = factor(number_recovered$variable, 
                                   levels = c("Only_from_gc","Only_from_imputation", "Two_consistent", "Two_inconsistent"),
                                   labels = c("only from genotype caller", "only from imputation", "called consistently from the two", "called inconsistently from the two"))
g_recovered = ggplot() + 
  geom_bar(data = number_recovered, aes(x = sample, y=value, fill = variable), stat = 'identity')  + 
  ylab('Fraction of variants from genotype data') + 
  xlab('Samples') + 
  scale_fill_brewer(palette = 'Set1') + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE, title=""))
ggsave(paste0("results/Number_recovered_variants_with_Imputation_",pattern,".png"), g_recovered, height = 4.2, width = 6.5)








### variants
readin_variants_df <- function(pattern){
  files = list.files('data/performance/')
  files = files[grepl("precision_recall_r", files)]
  files = files[grepl(pattern, files)]
  dat = NULL
  for(f in files){
    dati = read.table(paste0('data/performance/', f), sep=',', header = T, stringsAsFactors = F)
    sample = strsplit(f, '_')[[1]][1]
    dati$sample = sample
    if(grepl('inconsistent_prediction_number_precision_recall_r_', f)){
      tp = strsplit(f, '_')[[1]]
      mdl = paste(tp[length(tp)-1], tp[length(tp)],sep = '_')
      mdl = gsub('.txt', '', mdl)
      dati$Model = mdl
    }
    if(!'Model' %in% colnames(dati)){
      dati$Model = 'None'
    }
    if(!'Distance_to_peaks' %in% colnames(dati)){
      dati$Distance_to_peaks = 'no_threshold'
    }
    dat = rbind(dat, dati)
  }
  dat$pattern = pattern
  return(dat)
}

dat = NULL
for(minDP in c(2,3, 4,6)){
  dat = rbind(dat, readin_variants_df(paste0('minDP', minDP)))
}

dat$Group = paste0(dat$Group,'_', dat$Model)
dat$Group = gsub("_None", "", dat$Group)



### Genotype caller and Imputation
colList =  c("Precision", "Recall", "Called", "Pearson_correlation", "True_GT", "pattern")
metrics = dat[, c("sample", colList, "Group", "Distance_to_peaks")]
metrics = melt(metrics, id.vars = c("sample", "Group", "True_GT", "pattern", "Distance_to_peaks"))

metrics = metrics[metrics$Group %in% c("genotype_caller", "imputation"), ]
metrics$Group = factor(metrics$Group,  levels = c("genotype_caller", "imputation"),
                       labels = c("Genotype caller", "Imputation"))
metrics$variable = factor(metrics$variable,  levels = c("Called", "Recall", "Precision", "Pearson_correlation"),
                       labels = c("Number of called variants", "Recall", "Precision", "Pearson correlation with true dosage"))
metrics$True_GT = factor(metrics$True_GT, levels = c(0, 1, 2), labels = c("AA", "AB", "BB"))
metrics$Distance_to_peaks = factor(metrics$Distance_to_peaks)

g_metrics = ggplot(data = metrics[metrics$variable == 'Pearson correlation with true dosage', ], aes(x = pattern, y=value, fill = Group, color = Group)) + 
  facet_wrap(~ Distance_to_peaks, ncol = 3)  + 
  geom_boxplot()+ 
  xlab('Category of called variants') +
  ylab('Pearson correlation bewteen dosage') + 
  scale_fill_brewer(palette = 'Set3') + 
  scale_color_brewer(palette = 'Set3') + 
  guides(fill=guide_legend(title="")) + 
  theme_bw() + 
  guides(color = FALSE)
ggsave(paste0("results/Genotype_metrics_with_Imputation_Pearson.png"), g_metrics, height = 3, width = 5)


## Recall & Precision by minDP & distance to peaks
## Update (1/19/2021)
plot_precision_recall_gc_imputation <- function(metric, dat_plot, save_text){
  df = dat_plot[dat_plot$variable == metric, ]
  df$pattern = gsub('minDP', '', df$pattern)
  df$pattern = as.numeric(df$pattern)
  
  median = df %>% group_by(pattern, True_GT, Distance_to_peaks, Group) %>% summarise(median = median(value))
  sde = df %>% group_by(pattern, True_GT, Distance_to_peaks, Group) %>% summarise(sde = sd(value))
  median_sed = merge(median, sde, by = c("pattern", "True_GT", "Distance_to_peaks", "Group"))
  
  g_metrics = ggplot(data = median_sed, aes(x = pattern, y=median, color = Distance_to_peaks)) + 
    geom_point()+
    geom_line()+
    facet_wrap(~  Group + True_GT, ncol = 3)  + 
    xlab('Category of called variants (minDP threshold)') +
    ylab(metric) + 
    scale_color_viridis(discrete = TRUE) + 
    scale_fill_viridis(discrete = TRUE)+
    guides(fill=guide_legend(title="Distance to the nearest peaks")) + 
    theme_bw() + 
    geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = Distance_to_peaks), alpha = 0.3, color = "NA") + 
    guides(fill = FALSE)
  
  ggsave(paste0("results/",metric,"_",save_text, ".png"), g_metrics, height = 5.5, width = 10)
}


dat_plot = metrics[metrics$variable == 'Number of called variants', ]
plot_precision_recall_gc_imputation("Number of called variants", dat_plot, 'with_Imputation')

dat_plot = metrics[(metrics$variable != 'Number of called variants') & (metrics$variable != 'Pearson correlation with true dosage'), ]
plot_precision_recall_gc_imputation("Recall", dat_plot, 'with_Imputation')
plot_precision_recall_gc_imputation("Precision", dat_plot, 'with_Imputation')


## Inconsistent performance
colList =  c("Precision", "Recall", "Called", "Pearson_correlation", "True_GT")
metrics = dat[, c("sample", colList, "Group", "pattern", "Model")]
metrics = melt(metrics, id.vars = c("sample", "Group", "True_GT", "pattern", "Model"))

metrics = metrics[metrics$Group %in% c("inconsistent_called", "inconsistent_imputation", 
                                       "inconsistent_predicted_linear_regression", "inconsistent_predicted_logistic_regression", 
                                       "inconsistent_predicted_random_forest"), ]
metrics$Group = factor(metrics$Group,  
                       levels = c("inconsistent_called", "inconsistent_imputation",  
                                  "inconsistent_predicted_linear_regression", "inconsistent_predicted_logistic_regression", "inconsistent_predicted_random_forest"),
                       labels = c("Use Genotype caller", "Use Imputation", 
                                  "Predictions (linear regression)", 
                                  "Predictions (logistic regression)", 
                                  "Predictions (random forest)"))

metrics$variable = factor(metrics$variable,  levels = c("Called", "Recall", "Precision", "Pearson_correlation"),
                          labels = c("Number of called variants", "Recall", "Precision", "Pearson correlation with true dosage"))
metrics$True_GT = factor(metrics$True_GT, levels = c(0, 1, 2), labels = c("AA", "AB", "BB"))

plot_precision_recall_inconsistent <- function(metric, dat_plot, save_text){
  df = dat_plot[dat_plot$variable == metric, ]
  df$pattern = gsub('minDP', '', df$pattern)
  df$pattern = as.numeric(df$pattern)
  
  median = df %>% group_by(pattern, True_GT, Group) %>% summarise(median = median(value))
  sde = df %>% group_by(pattern, True_GT, Group) %>% summarise(sde = sd(value))
  median_sed = merge(median, sde, by = c("pattern", "True_GT", "Group"))
  
  g_metrics = ggplot(data = median_sed, aes(x = pattern, y=median, color = Group)) + 
    geom_point()+
    geom_line()+
    facet_wrap(~  True_GT, ncol = 3)  + 
    xlab('Category of called variants (minDP threshold)') +
    ylab(metric) + 
    scale_color_viridis(discrete = TRUE) + 
    scale_fill_viridis(discrete = TRUE)+
    guides(fill=guide_legend(title="Distance to the nearest peaks")) + 
    theme_bw() + 
    geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = Group), alpha = 0.3, color = "NA") + 
    guides(fill = FALSE)
  
  ggsave(paste0("results/",metric,"_",save_text, ".png"), g_metrics, height = 2.3, width = 10)
}

dat_plot = metrics[metrics$variable == 'Number of called variants', ]
plot_precision_recall_inconsistent("Number of called variants", dat_plot, save_text = 'Inconsistent')

dat_plot = metrics[metrics$variable == 'Pearson correlation with true dosage', ]
plot_precision_recall_inconsistent("Pearson correlation with true dosage", dat_plot, save_text = 'Inconsistent')


dat_plot = metrics[(metrics$variable != 'Number of called variants') & (metrics$variable != 'Pearson correlation with true dosage'), ]
plot_precision_recall_inconsistent("Recall", dat_plot, save_text = 'Inconsistent')
plot_precision_recall_inconsistent("Precision", dat_plot, save_text = 'Inconsistent')





# Combined performance
colList =  c("Precision", "Recall", "Called", "Pearson_correlation", "True_GT", "pattern")
dat_all = dat[dat$Distance_to_peaks == 'no_threshold', ]
metrics = dat_all[, c("sample", colList, "Group")]
metrics = melt(metrics, id.vars = c("sample", "Group", "True_GT", "pattern"))

metrics = metrics[metrics$Group %in% c("combined_no_inconsistent", "combined_with_inconsistent_random_forest", 
                                       "combined_with_inconsistent_linear_regression", "combined_with_inconsistent_logistic_regression",
                                       "imputation", "genotype_caller"), ]
metrics$Group = factor(metrics$Group,  
                       levels = c("genotype_caller", "imputation", "combined_no_inconsistent", 
                                  "combined_with_inconsistent_random_forest", "combined_with_inconsistent_linear_regression", "combined_with_inconsistent_logistic_regression"),
                       labels = c("Genotype caller", "Imptuation","Combined (remove inconsistent calls)", 
                                  "Combined (predictions from random forest)", "Combined (predictions from linear regression)", "Combined (predictions from logistic regression)"))

metrics$variable = factor(metrics$variable,  levels = c("Called", "Recall", "Precision", "Pearson_correlation"),
                          labels = c("Number of called variants", "Recall", "Precision", "Pearson correlation with true dosage"))
metrics$True_GT = factor(metrics$True_GT, levels = c(0, 1, 2), labels = c("AA", "AB", "BB"))


dat_plot = metrics[metrics$variable == 'Number of called variants', ]
plot_precision_recall_inconsistent("Number of called variants", dat_plot, save_text = 'combined')

dat_plot = metrics[metrics$variable == 'Pearson correlation with true dosage', ]
plot_precision_recall_inconsistent("Pearson correlation with true dosage", dat_plot, save_text = 'combined')


dat_plot = metrics[(metrics$variable != 'Number of called variants') & (metrics$variable != 'Pearson correlation with true dosage'), ]
plot_precision_recall_inconsistent("Recall", dat_plot, save_text = 'combined')
plot_precision_recall_inconsistent("Precision", dat_plot, save_text = 'combined')





## Integrate over all genotypes
library(dplyr)
dat_summary = dat_all %>% group_by(sample, Group, Distance_to_peaks, pattern) %>% summarise(Golden_standard = sum(Golden_standard), 
                                                                                        Called = sum(Called), 
                                                                                        Called_and_True = sum(Called_and_True),
                                                                                        Pearson_correlation = mean(Pearson_correlation))
dat_summary$Precision = dat_summary$Called_and_True / dat_summary$Called
dat_summary$Recall = dat_summary$Called_and_True / dat_summary$Golden_standard
dat_summary = dat_summary[dat_summary$Group %in% c("combined_no_inconsistent", "combined_with_inconsistent_random_forest", 
                                       "combined_with_inconsistent_linear_regression", "combined_with_inconsistent_logistic_regression",
                                       "imputation", "genotype_caller"), ]

# add true numbers
#dat_true = dat_summary
#dat_true$Called = dat_true$Golden_standard
#dat_true$Called_and_True = dat_true$Golden_standard
#dat_true$Group = "True"
#dat_true$Precision = 1
#dat_true$Recall = 1
#dat_true$Pearson_correlation = 1
#dat_true = dat_true[!duplicated(dat_true), ]

#dat_summary = as.data.frame(rbind(dat_summary, dat_true))

dat_summary$Group = factor(dat_summary$Group,  
                           levels = c("True","genotype_caller", "imputation", "combined_no_inconsistent", 
                                      "combined_with_inconsistent_random_forest", "combined_with_inconsistent_linear_regression", "combined_with_inconsistent_logistic_regression"),
                           labels = c("True","Genotype caller", "Imptuation","Combined (remove inconsistent calls)", 
                                      "Combined (predictions from random forest)", "Combined (predictions from linear regression)", "Combined (predictions from logistic regression)"))

plot_precision_recall_overall <- function(metric, df, save_text){
  df[, "value"] = df[, metric]
  df$pattern = gsub('minDP', '', df$pattern)
  df$pattern = as.numeric(df$pattern)
  
  median = df %>% group_by(pattern, Group) %>% summarise(median = median(value))
  sde = df %>% group_by(pattern, Group) %>% summarise(sde = sd(value))
  median_sed = merge(median, sde, by = c("pattern", "Group"))
  
  g_metrics = ggplot(data = median_sed, aes(x = pattern, y=median, color = Group)) + 
    geom_point()+
    geom_line()+
    xlab('Category of called variants (minDP threshold)') +
    ylab(metric) + 
    scale_color_viridis(discrete = TRUE) + 
    scale_fill_viridis(discrete = TRUE)+
    guides(fill=guide_legend(title="Distance to the nearest peaks")) + 
    theme_bw() + 
    geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = Group), alpha = 0.3, color = "NA") + 
    guides(fill = FALSE)
  
  ggsave(paste0("results/",metric,"_",save_text, ".png"), g_metrics, height = 2.3, width = 6)
}


plot_precision_recall_overall("Recall", dat_summary, "combined_allGT")
plot_precision_recall_overall("Precision", dat_summary, "combined_allGT")
plot_precision_recall_overall("Pearson_correlation", dat_summary, "combined_allGT")

