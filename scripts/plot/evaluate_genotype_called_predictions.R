library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(cowplot)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

training_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model'
dat_all = NULL
dat_GT = NULL
dat_distance = NULL

for(SUFFIX in c("_THR0.0.txt","_THR0.1.txt", "_THR0.2.txt")){
  for(minDP in c(2, 3, 4, 5, 6,8,10)){
    saved_file = paste0('data/Model/minDP', minDP, SUFFIX, '.RData')
    if( file.exists(saved_file)){
      load(saved_file)
      resulti[[1]]$THR = gsub("_", "",gsub(".txt", "", SUFFIX))
      resulti[[2]]$THR = gsub("_", "",gsub(".txt", "", SUFFIX))
      resulti[[3]]$THR = gsub("_", "",gsub(".txt", "", SUFFIX))
      dat_all = rbind(dat_all, resulti[[1]])
      dat_GT  = rbind(dat_GT, resulti[[2]])
      dat_distance = rbind(dat_distance, resulti[[3]])
    }else{
      print(paste0(saved_file, ' does not exist'))
    }
  }
}

dat_all = dat_all[, colnames(dat_all)[colnames(dat_all) != 'OR']]


### Plot
metrics = melt(dat_all[dat_all$Metric != 'Pearson', ], id.vars = c("sample", "pattern", "Metric", "Group", "THR"))

metrics$value = as.numeric(metrics$value)
metrics$variable = factor(metrics$variable,  levels = c("GC", "IP",  "LiR", "LgR", "OR_R", 'RF'),
                          labels = c("Genotype caller","Imputation", 
                                     "Linear regression pred.", "Logistic regression pred.", "Ordinal regression pred.", "Random forest pred."))

median = metrics %>% group_by(pattern, variable, Metric, THR) %>% summarise(median = median(value))
sde = metrics %>% group_by(pattern, variable, Metric, THR) %>% summarise(sde = sd(value))
median_sed = merge(median, sde, by = c("pattern", "variable", "Metric", "THR"))
median_sed$pattern = as.numeric(gsub('minDP', '', median_sed$pattern))




# Only GC and Imputation
dat_plot = median_sed[(median_sed$variable == "Genotype caller") | (median_sed$variable == "Imputation"), ]

g_metrics = ggplot(data = dat_plot[dat_plot$Metric == 'MSE', ], aes(x = pattern, y=median, color = variable)) +
  #facet_wrap(~THR) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('MSE') +
  scale_color_manual(values = c("red", "blue")) +
  #scale_color_viridis(discrete = TRUE) +
  #scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = variable), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_MSE_part.png"), g_metrics, height = 2.5, width = 5)



g_metrics = ggplot(data = dat_plot[dat_plot$Metric == 'Spearman', ], aes(x = pattern, y=median, color = variable)) +
  #facet_wrap(~THR) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('Spearman correlation') +
  scale_color_manual(values = c("red", "blue")) +
  #scale_color_viridis(discrete = TRUE) +
  #scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = variable), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_Spearman_part.png"), g_metrics, height = 2.5, width = 5)



# All
g_metrics = ggplot(data = median_sed[median_sed$Metric == 'MSE', ], aes(x = pattern, y=median, color = variable)) +
  #facet_wrap(~THR) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('MSE') +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = variable), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_MSE.png"), g_metrics, height = 2.5, width = 5)





g_metrics = ggplot(data = median_sed[median_sed$Metric == 'Spearman', ], aes(x = pattern, y=median, color = variable)) +
  #facet_wrap(~THR) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('Spearman correlation') +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = variable), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_Spearman.png"), g_metrics, height = 2.5, width = 5)




### Plot out the linear regression

coef = read.table('data/Model/linear_regression_coef.txt', sep='\t', header = T)
coef = melt(coef, id.vars = 'Read_depth_filter')
g = ggplot(data = coef) + 
  geom_tile(aes(x = Read_depth_filter, y = variable, fill=value)) + 
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab" ) + 
  ylab('Coefficient')

ggsave(paste0("results/Linear_regression_coef.pdf"), g, height = 3.5, width = 5)


## Number of variants
dat_example = read.table('data/Performance/Overall_performance_HG00105.txt', sep='\t', header = T)
dat_example = dat_example[, c("minDP", "Number", "method")]
dat_example$method = factor(dat_example$method, levels = c("GC", "Imputation", "Imputed_dosage", "Dosage", 
                                                         "Y_linear_regression", "Y_logistic_regression", "Y_ordinal_regression", "Y_random_forest"),
                           labels = c("Only genotype caller", "Only imputation", "Both", "Both",
                                      "Both", "Both", "Both", "Both"))
dat_example = dat_example[!duplicated(dat_example), ]
dat_example$minDP = as.factor(dat_example$minDP)
g = ggplot(data = dat_example) + geom_bar(aes(x = minDP, y = Number, fill = method), stat = 'identity') +
  ylab('Number of called variants')
ggsave(paste0("results/dat_example.pdf"), g, height = 3.5, width = 5)

