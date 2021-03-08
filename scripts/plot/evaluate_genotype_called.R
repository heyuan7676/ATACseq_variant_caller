library("ggsci")
setwd('/Users/Yuan/Desktop/Variant_calling/')

library(ggplot2)
library(reshape2)

dat_all = NULL
dat_SNPs = NULL
for(minDP in seq(2, 9)){
  token_1 = read.table(paste0('data/performance/merge_left/All_variants_minDP', minDP, '.txt'), sep='\t', header=T)
  dat_all = rbind(dat_all, token_1)
  
  token_2 = read.table(paste0('data/performance/merge_left/All_SNPs_minDP', minDP,'.txt'), sep='\t', header=T)
  dat_SNPs = rbind(dat_SNPs, token_2)
}

dat_all$pool = 'All variants'
dat_SNPs$pool = 'All SNPs' 
dat = rbind(dat_all, dat_SNPs)
dat$minDP = factor(dat$minDP)

### Number of recovered variants
number_recovered = dat[, c("minDP", "N_variants_by_ATAC", "pool")]
number_recovered = melt(number_recovered, id.vars = c("minDP", "pool"))
g_recovered = ggplot() + 
  geom_boxplot(data = number_recovered, aes(x = minDP, y=value, fill = pool), lwd=0.2, outlier.size = 0.2)  + 
  ylab('Number of variants called from ATAC-seq reads') + 
  xlab('minimum read depth for each variant across samples') + 
  scale_fill_brewer(palette = 'Set1') + 
  guides(fill=guide_legend(title="")) + 
  theme_bw()
ggsave("results/Number_recovered_variants.png", g_recovered, height = 4, width = 8)



### Percentage of recovered variants
recovered = dat[,c("minDP","Recovered_percentage", "Correctly_performance_percentage", "pool")]
colnames(recovered) = c("minDP", "Called from ATAC-seq", "Correctly called from ATAC-seq", "pool")
recovered = melt(recovered, id.vars = c("minDP", "pool"))
#recovered$minDP = factor(recovered$minDP)

g_percentage_recovered = ggplot() + 
  geom_boxplot(data = recovered, aes(x = minDP, y=value, fill = pool), lwd=0.2, outlier.size = 0.2)  + 
  ylab('Fraction of variants from genotype data') + 
  xlab('minimum read depth for each variant across samples') + 
  scale_fill_brewer(palette = 'Greens') + 
  guides(fill=guide_legend(title="")) +
  facet_wrap(~variable, ncol = 2) + 
  theme_bw()
ggsave("results/Factorion_recovered_variants.png", g_percentage_recovered, height = 3, width = 10)



### Among the tested variants:

recall = dat[,c("minDP","Sens_AA", "Sens_AB", "Sens_BB", "pool")]
colnames(recall) = c("minDP", "AA", "AB", "BB", "pool")
recall = melt(recall, id.vars = c("minDP", "pool"))
#recovered$minDP = factor(recovered$minDP)

g_recall = ggplot(data = recall, aes(x = minDP, y=value, fill = variable)) + 
  geom_bar(stat='summary', fun.data = mean_se)  + 
  geom_errorbar(stat='summary', width=.2, fun.data = mean_se) + 
  ylab('Recall of the tested variants') + 
  xlab('minimum read depth for each variant across samples') + 
  scale_fill_brewer(palette = 'Blues') + 
  guides(fill=guide_legend(title="Genotype")) +
  facet_wrap(~pool + variable, nrow=2) +
  theme_bw()
ggsave("results/Genotype_recall.png", g_recall)



precision = dat[,c("minDP","Spec_AA", "Spec_AB", "Spec_BB", "pool")]
colnames(precision) = c("minDP", "AA", "AB", "BB", "pool")
precision = melt(precision, id.vars = c("minDP", "pool"))
#recovered$minDP = factor(recovered$minDP)

g_precision = ggplot(data = precision, aes(x = minDP, y=value, fill = variable)) + 
  geom_bar(stat='summary', fun.data = mean_se)  + 
  geom_errorbar(stat='summary', width=.2, fun.data = mean_se) + 
  ylab('Precision of the tested variants') + 
  xlab('minimum read depth for each variant across samples') + 
  scale_fill_brewer(palette = 'Dark2') + 
  guides(fill=guide_legend(title="Genotype")) +
  facet_wrap(~pool + variable, nrow = 2) + 
  theme_bw()

ggsave("results/Genotype_precision.png", g_precision)


## compare SNPs with all variants
cols = c("Sample", "minDP", "N_variants_by_ATAC", "Correctly_performance_percentage")
dat_compare = merge(dat_SNPs[, cols], dat_all[, cols], by = c("Sample", "minDP"))

min(dat_compare$N_variants_by_ATAC.x / dat_compare$N_variants_by_ATAC.y)
max(dat_compare$N_variants_by_ATAC.x / dat_compare$N_variants_by_ATAC.y)

dat_compare$Correctly_performance_percentage.x / dat_compare$Correctly_performance_percentage.y

t.test(dat_compare$Correctly_performance_percentage.x, dat_compare$Correctly_performance_percentage.y, paired = TRUE)
min((dat_compare$Correctly_performance_percentage.x - dat_compare$Correctly_performance_percentage.y) / dat_compare$Correctly_performance_percentage.y)
max((dat_compare$Correctly_performance_percentage.x - dat_compare$Correctly_performance_percentage.y) / dat_compare$Correctly_performance_percentage.y)








### imputed vs. original called
files = list.files('data/performance/combined/')
files = files[grepl("called_in_both", files)]
files = files[grepl("allVariants", files)]
dat = NULL
for(f in files){
  dati = read.table(paste0('data/performance/combined/', f), sep='\t', header = T, stringsAsFactors = F)
  dat = rbind(dat, dati)
}


### Number of recovered variants
for(col in c("only_in_orignal", "called_in_both_consistent", "called_in_both_inconsistently", "only_in_imputed")){
  dat[,col] = dat[,col] / dat$Variants_in_WGS
}
number_recovered = dat[, c("sample", c("only_in_orignal", "called_in_both_consistent", "called_in_both_inconsistently", "only_in_imputed"))]
number_recovered = melt(number_recovered, id.vars = c("sample"))
number_recovered$variable = factor(number_recovered$variable, 
                                   levels = c("only_in_orignal","only_in_imputed", "called_in_both_consistent", "called_in_both_inconsistently"),
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
ggsave("results/Number_recovered_variants_with_Imputation.png", g_recovered, height = 4.2, width = 6.5)








### imputed vs. original called
files = list.files('data/performance/combined/')
files = files[grepl("performance_one_method", files)]
files = files[grepl("allVariants", files)]
dat = NULL
for(f in files){
  dati = read.table(paste0('data/performance/combined/', f), sep='\t', header = T, stringsAsFactors = F)
  dat = rbind(dat, dati)
}


### Recall and Precision
colList =  c("Recall_AA_all", "Recall_AB_all", "Recall_BB_all", "Precision_AA", "Precision_AB", "Precision_BB")
metrics = dat[, c("Sample", colList, "group")]
metrics = melt(metrics, id.vars = c("Sample", "group"))
metrics$metric = sapply(metrics$variable, function(x) strsplit(as.character(x), '_')[[1]][1])
metrics$GT = sapply(metrics$variable, function(x) strsplit(as.character(x), '_')[[1]][2])

metrics = metrics[metrics$group %in% c("performance_only_in_original", "performance_consistent", "performance_only_in_imputed"), ]
metrics$group = factor(metrics$group,  levels = c("performance_only_in_original", "performance_only_in_imputed", "performance_consistent"),
                       labels = c("only from Genotype caller", "only from Imputation", "consistent call from both"))

g_metrics = ggplot(data = metrics, aes(x = group, y=value, fill = group)) + 
  geom_bar(stat='summary', fun.data = mean_se)  + 
  geom_errorbar(stat='summary', width=.2, fun.data = mean_se) + 
  facet_wrap(~metric + GT)  + 
  xlab('Category of called variants') +
  ylab('Performance metric') + 
  scale_fill_brewer(palette = 'Set1') + 
  guides(fill=guide_legend(title="")) + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) 
ggsave("results/Genotype_metrics_with_Imputation.png", g_metrics, height = 5, width = 10)



## Overall performance
colList =  c("Precision_overall", "Recall_overall")
metrics = dat[, c("Sample", colList, "group", "minDP")]
metrics = melt(metrics, id.vars = c("Sample", "group", "minDP"))
metrics$variable = factor(metrics$variable,  levels = colList)

metrics = metrics[metrics$group %in% c("performance_in_original", "performance_in_imputed", 
                                       "performance_overall_noInconsistent","performance_overall_withPP"), ]
metrics$group = factor(metrics$group,  
                       levels = c("performance_in_original", "performance_in_imputed", 
                                  "performance_overall_noInconsistent", "performance_overall_withPP"),
                       labels = c("Genotype caller", "Imputation", "Combined without inconsistent calls", "Combined with inconsistent calls"))

metrics$minDP = factor(metrics$minDP)

g_metrics = ggplot(data = metrics, aes(x = group, y=value, fill = group)) + 
  geom_bar(stat='summary', fun.data = mean_se)  + 
  geom_errorbar(stat='summary', width=.2, fun.data = mean_se) + 
  facet_wrap(~variable + minDP)  + 
  ylab('Metrics') + 
  xlab('Category of called variants') + 
  scale_fill_brewer(palette = 'Paired') + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE, title=""))
ggsave("results/Genotype_performance_with_Imputation_withPP.png", g_metrics, height = 5, width = 6)











## for the in-consistent part
files = list.files('data/performance/combined/')
files = files[grepl("called_inconsistently", files)]
files = files[grepl("allVariants", files)]
dat_incons = NULL
for(f in files){
  dati = read.table(paste0('data/performance/combined/', f), sep='\t', header = T, stringsAsFactors = F)
  dat_incons = rbind(dat_incons, dati)
}

dat_incons$PP = round(dat_incons$PP, 2)

dat_incons$AA_BB_to_AB_recall = dat_incons$AA_BB_to_AB / dat_incons$Number_AB_all
dat_incons$AB_to_AA_BB_recall = dat_incons$AB_to_AA_BB / (dat_incons$Number_AA_all + dat_incons$Number_BB_all)
dat_incons$AA_to_BB_recall = dat_incons$AA_to_BB / (dat_incons$Number_AA_all + dat_incons$Number_BB_all)

dat_incons_recall_cols = colnames(dat_incons)[grep('recall', colnames(dat_incons))]
dat_incons_recall = dat_incons[, c(c("sample"), dat_incons_recall_cols, "PP")]
dat_incons_recall = melt(dat_incons_recall, id.vars = c("sample", "PP"))
dat_incons_recall$X = 1
dat_incons_recall = as.data.frame(rbind(dat_incons_recall, c("HG00096", "0.33", "AA_BB_to_AB_recall", 0, 0)))
dat_incons_recall = as.data.frame(rbind(dat_incons_recall, c("HG00096", "0.33", "AA_BB_to_AB_recall", 0, 2)))
dat_incons_recall$variable = factor(dat_incons_recall$variable, 
                                    levels = c("AA_BB_to_AB_recall", "AB_to_AA_BB_recall", "AA_to_BB_recall"),
                                    labels = c("AA/BB -> AB", "AB -> AA/BB", "AA/BB -> BB/AA"))
dat_incons_recall$PP = factor(dat_incons_recall$PP, 
                              levels = c("0.33", "0.5", "1"), 
                              labels = c("PP=0.33", "PP=0.5", "PP=1"))
dat_incons_recall$value = as.numeric(dat_incons_recall$value)

metrics_incons = ggplot(data = dat_incons_recall, aes(x = X, y=value)) + 
  facet_wrap(~variable + PP) + 
  geom_bar(stat='summary', fun.data = mean_se, position = 'dodge',  fill = "#999999", width = 0.8)  + 
  geom_errorbar(stat='summary', width=0.1, fun.data = mean_se, position = 'dodge') + 
  ylab('Recall') + 
  xlab('Posterior probability from genotype caller') + 
  guides(fill=guide_legend(title="")) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("results/Genotype_metrics_inconsistent_recall.png", metrics_incons, height = 4, width = 4)




dat_incons_precision_cols = colnames(dat_incons)[grep('precision', colnames(dat_incons))]
dat_incons_precision = dat_incons[, c(c("sample"), dat_incons_precision_cols, "PP")]
dat_incons_precision = melt(dat_incons_precision, id.vars = c("sample", "PP"))
dat_incons_precision[!complete.cases(dat_incons_precision$value),"value"] = 0

dat_incons_precision$group = "group"
dat_incons_precision[grep("AA_BB_to_AB", dat_incons_precision$variable), "group"] = "AA_BB_to_AB"
dat_incons_precision[grep("AB_to_AA_BB", dat_incons_precision$variable), "group"] = "AB_to_AA_BB"
dat_incons_precision[grep("AA_to_BB", dat_incons_precision$variable), "group"] = "AA_to_BB"

dat_incons_precision$useValue = "value"
dat_incons_precision[grep("original", dat_incons_precision$variable), "useValue"] = "genotype caller"
dat_incons_precision[grep("imputed", dat_incons_precision$variable), "useValue"] = "imputed"

dat_incons_precision$group = factor(dat_incons_precision$group, 
                                    levels = c("AA_BB_to_AB", "AB_to_AA_BB", "AA_to_BB"),
                                    labels = c("AA/BB -> AB","AB -> AA/BB","AA/BB -> BB/AA"))
dat_incons_precision$useValue = factor(dat_incons_precision$useValue)

dat_incons_precision$PP = factor(dat_incons_precision$PP, levels = c("0.33", "0.5", "1"), labels = c("PP=0.33", "PP=0.5", "PP=1"))

metrics_incons_precision = ggplot(data = dat_incons_precision, aes(x = useValue, y=value, fill = useValue)) + 
  facet_wrap(~ group + PP) + 
  geom_bar(stat='summary', fun.data = mean_se, width = 0.8)  + 
  geom_errorbar(stat='summary', width=.9, fun.data = mean_se) + 
  ylab('Precision') + 
  xlab('Category of the called variants') + 
  scale_fill_manual(values = c("red", "#56B4E9")) + 
  guides(fill=guide_legend(title="")) + 
  theme_bw() +
  theme(axis.text.x = element_blank())
ggsave("results/Genotype_metrics_inconsistent_precision.png", metrics_incons_precision, height = 5, width = 6.5)








dat = data.frame("group" = c("only_in_orignal", "called_in_both_consistent", "called_in_both_inconsistently", "only_in_imputed"),
             "value" = c(134893, 1341079, 380586, 3506902))
dat$group = factor(dat$group,
                   levels = rev(c("only_in_imputed", "called_in_both_consistent", "called_in_both_inconsistently", "only_in_orignal")))

g_plot = ggplot() + 
  geom_bar(aes(x = "HG00096", y=value, fill=group), 
           position = 'stack', 
           stat = 'identity', 
           data = dat)  + 
  xlab("") + 
  ylab("Number of SNPs") + 
  scale_fill_brewer(palette = 'Set3') + 
  coord_flip() + 
  theme_bw() + 
  theme(panel.border = element_blank())

ggsave('results/HG00096_compare_original_imputed.png', g_plot, width = 10, height = 2)

