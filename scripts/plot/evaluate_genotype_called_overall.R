library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(cowplot)

training_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model'
files = list.files('data/Performance')
files = files[grepl('Overall_performance', files)]
files = files[!grepl('byCategory', files)]
files = files[!grepl('Gencove', files)]

dat = NULL
for(f in files){
  print(f)
  dt = read.table(paste0('data/Performance/', f), sep='\t', header=T)
  dat = rbind(dat, dt)
}

dt = read.table('data/Performance/Overall_performance_Gencove.txt', sep='\t', header=T)
for(minDP in seq(2, 10)){
  dt$minDP = minDP
  dat = rbind(dat, dt)
}

dat_number = dat[,c("sample", "minDP", "Number", "method")]
dat_number$method = factor(dat_number$method, levels = c("GC", "Imputation", "Imputed_dosage", "Dosage", 
                                           "Y_linear_regression", "Y_logistic_regression", "Y_ordinal_regression", "Y_random_forest", "Gencove"),
                    labels = c("Genotype caller", "Imputation", "Combined", "Combined",
                               "Combined", "Combined", "Combined", "Combined", "Gencove"))

dat_number = dat_number[!duplicated(dat_number),]
median = dat_number %>% group_by(method, minDP) %>% summarise(median = median(Number))
sde = dat_number %>% group_by(method, minDP) %>% summarise(sde = sd(Number))

median_sed = merge(median, sde, by = c("method", "minDP"))

g_metrics = ggplot(data = median_sed, aes(x = minDP, y=median, color = method)) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('Number of called variants') +
  scale_color_manual(values = c("red", "blue", "green", "purple")) + 
  #scale_color_viridis(discrete = TRUE) +
  #scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  #geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = method), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_overall_Number.png"), g_metrics, height = 2.5, width = 5)





dat = dat[dat$method != 'GC', ]
dat = dat[dat$method != 'Imputation', ]
dat$method = factor(dat$method, levels = c("GC", "Imputation", "Imputed_dosage", "Dosage", 
                                           "Y_linear_regression", "Y_logistic_regression", "Y_ordinal_regression", "Y_random_forest",
                                           "Gencove"),
                    labels = c("Only genotype caller", "Only imputation", "Combined (use Imputation)", "Combined (use genotype caller)",
                               "Combined (linear regression pred.)", "Combined (logistic regression pred.)", 
                               "Combined (ordinal regression pred.)", "Combined (random forest pred.)", "Gencove"))


### Spearman correlation
median = dat %>% group_by(method, minDP) %>% summarise(median = median(correlation))
sde = dat %>% group_by(method, minDP) %>% summarise(sde = sd(correlation))

median_sed = merge(median, sde, by = c("method", "minDP"))

g_metrics = ggplot(data = median_sed, aes(x = minDP, y=median, color = method)) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('Spearman correlation') +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  #geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = method), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_overall_spearman_correlation.png"), g_metrics, height = 2.5, width = 5.5)




### MSE
median = dat %>% group_by(method, minDP) %>% summarise(median = median(MSE))
sde = dat %>% group_by(method, minDP) %>% summarise(sde = sd(MSE))

median_sed = merge(median, sde, by = c("method", "minDP"))

g_metrics = ggplot(data = median_sed, aes(x = minDP, y=median, color = method)) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('MSE') +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  #geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = method), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_overall_MSE.png"), g_metrics, height = 2.5, width = 5.5)




### MSE - AA
median = dat %>% group_by(method, minDP) %>% summarise(median = median(MSE_AA))
sde = dat %>% group_by(method, minDP) %>% summarise(sde = sd(MSE_AA))

median_sed = merge(median, sde, by = c("method", "minDP"))

g_metrics = ggplot(data = median_sed, aes(x = minDP, y=median, color = method)) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('MSE - AA') +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  #geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = method), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_overall_MSE_AA.png"), g_metrics, height = 2.5, width = 5)




### MSE - AB
median = dat %>% group_by(method, minDP) %>% summarise(median = median(MSE_AB))
sde = dat %>% group_by(method, minDP) %>% summarise(sde = sd(MSE_AB))

median_sed = merge(median, sde, by = c("method", "minDP"))

g_metrics = ggplot(data = median_sed, aes(x = minDP, y=median, color = method)) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('MSE - AB') +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  #geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = method), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_overall_MSE_AB.png"), g_metrics, height = 2.5, width = 5)




### MSE - BB
median = dat %>% group_by(method, minDP) %>% summarise(median = median(MSE_BB))
sde = dat %>% group_by(method, minDP) %>% summarise(sde = sd(MSE_BB))

median_sed = merge(median, sde, by = c("method", "minDP"))

g_metrics = ggplot(data = median_sed, aes(x = minDP, y=median, color = method)) +
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('MSE - BB') +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  #geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = method), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_overall_MSE_BB.png"), g_metrics, height = 2.5, width = 5)



### By category
training_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model'
files = list.files('data/Performance/')
files = files[grep('Overall_performance', files)]
files = files[grepl('_byCategory.txt', files)]
files = files[!grepl('minDP', files)]

dat = NULL
for(f in files){
  print(f)
  dt = read.table(paste0('data/Performance/', f), sep='\t', header=T)
  dat = rbind(dat, dt)
}



dat$Variant_category = factor(dat$Variant_category, levels = c("Dosage", "Imputed_dosage", 
                                           "Y_linear_regression", "Y_logistic_regression", "Y_ordinal_regression_R", "Y_random_forest"),
                    labels = c("Genotype caller","Imputation", 
                               "Linear regression pred.", "Logistic regression pred.", "Ordinal regression pred.", "Random forest pred."))


### Spearman correlation
median = dat %>% group_by(method, minDP, Variant_category) %>% summarise(median = median(correlation))
sde = dat %>% group_by(method, minDP, Variant_category) %>% summarise(sde = sd(correlation))

median_sed = merge(median, sde, by = c("method", "minDP", "Variant_category"))

g_metrics = ggplot(data = median_sed, aes(x = minDP, y=median, color = Variant_category)) +
  facet_wrap(~method) + 
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('Spearman correlation') +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  #geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = method), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_overall_spearman_correlation_byCat.png"), g_metrics, height = 3, width = 5.5)



### Spearman correlation
median = dat %>% group_by(method, minDP, Variant_category) %>% summarise(median = median(MSE))
sde = dat %>% group_by(method, minDP, Variant_category) %>% summarise(sde = sd(MSE))

median_sed = merge(median, sde, by = c("method", "minDP", "Variant_category"))

g_metrics = ggplot(data = median_sed, aes(x = minDP, y=median, color = Variant_category)) +
  facet_wrap(~method) + 
  geom_point()+
  geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('MSE') +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)+
  guides(fill=guide_legend(title="Method")) +
  theme_bw() +
  #geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = method), alpha = 0.3, color = "NA") +
  guides(fill = FALSE)

ggsave(paste0("results/Evaluation_metrics_overall_MSE_byCat.png"), g_metrics, height = 3, width = 5.5)

# Number
files = list.files('data/Performance/')
files = files[grepl('Overall_variants_number', files)]

dat_n = NULL
for(f in files){
  print(f)
  dt = read.table(paste0('data/Performance/', f), sep='\t', header=T)
  dat_n = rbind(dat_n, dt)
}

for(col in c("GC_only", "Imputation_only", "Both")){
  dat_n[,col] = dat_n[,col] / dat_n[, "WGS"]
}

dat_n$minDP = factor(dat_n$minDP)
dat_n = melt(dat_n, id.vars = c("sample", "minDP"))
dat_n = dat_n[dat_n$variable!='WGS', ]

dat_n$variable = factor(dat_n$variable, levels = c("GC_only", "Imputation_only", "Both"),
                        labels = c("Only from Genotype caller", 
                                   "Only from Imputation", 
                                   "Called from the two methods"))

#median = dat_n %>% group_by(minDP, variable) %>% summarise(median = median(value))
#sde = dat_n %>% group_by(minDP, variable) %>% summarise(sde = sd(value))

#median_sed = merge(median, sde, by = c( "minDP", "variable"))
#median_sed = median_sed[median_sed$variable != 'WGS', ]

g = ggplot(data = dat_n[dat_n$sample == 'HG00105', ], aes(x = minDP, y=value, fill = variable)) +
  geom_bar(stat = 'identity')+
  #geom_line()+
  xlab('Category of called variants (minDP threshold)') +
  ylab('Fraction of variants from WGS\n that are recovered') + 
  scale_fill_manual(values = c("red", "blue", "green")) +
  #scale_color_viridis(discrete = TRUE) +
  #scale_fill_viridis(discrete = TRUE)+
  #guides(fill=guide_legend(title="Method")) +
  theme_bw()
  #geom_ribbon(aes(ymax = median + sde, ymin = median - sde,  fill = method), alpha = 0.3, color = "NA") +
  #guides(fill = FALSE)

ggsave(paste0("results/dat_example.pdf"), g, height = 3.5, width = 5)

