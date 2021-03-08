setwd('/Users/Yuan/Desktop/Variant_calling/')

library(ggplot2)
library(reshape2)

dat_wgs = read.table('data/Evaluation_metrics/fastQTL_results_summary.txt', sep='\t', header = T, stringsAsFactors = F)
dat_wgs$cisDist = factor(dat_wgs$cisDist)
dat_wgs$peak_calling = factor(dat_wgs$peak_calling, levels = c("MACS2", "MACS2/combined", "Genrich", "Genrich/combined"),
                              labels = c("MACS2 (call-merge)", "MACS2 (merge-call)", "Genrich (call-merge)", "Genrich (merge-call)"))

g = ggplot(data = dat_wgs) + 
  geom_bar(aes(x=cisDist, y=QTLs, fill=peak_calling), stat = 'identity', position = "dodge") + 
  ylab('Number of QTLs at FDR < 0.05') + 
  xlab('Window size around Peak middle point to include variants') + 
  ggtitle('Results from fastQTL using true genotype') + 
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
  scale_fill_brewer(palette = "Paired")

ggsave(paste0('results/fastQTL_results_true_genotype.png'), height = 2.5, width = 5)



dat_wgs_2 = read.table('data/Evaluation_metrics/matrixEQTL_results_summary.txt', sep='\t', header = T, stringsAsFactors = F)
dat_wgs_2$cisDist = factor(dat_wgs_2$cisDist)
dat_wgs_2$peak_calling = factor(dat_wgs_2$peak_calling, levels = c("MACS2", "MACS2/combined", "Genrich", "Genrich/combined"),
                              labels = c("MACS2 (call-merge)", "MACS2 (merge-call)", "Genrich (call-merge)", "Genrich (merge-call)"))
g = ggplot(data = dat_wgs_2) + 
  geom_bar(aes(x=cisDist, y=QTLs, fill=peak_calling), stat = 'identity', position = "dodge") + 
  ylab('Number of QTLs at FDR < 0.05') + 
  xlab('Window size around Peak regions to include variants') + 
  ggtitle('Results from matrixEQTL using true genotype') + 
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
  scale_fill_brewer(palette = "Paired")

ggsave(paste0('results/matrixEQTL_results_true_genotype.png'), height = 2.5, width = 5)


colnames(dat_wgs_2) = colnames(dat_wgs)
dat_wgs$tool = 'fastQTL'
dat_wgs_2$tool = 'matrixEQTL'
dat_wgs_two = rbind(dat_wgs, dat_wgs_2)

dat_wgs_two = dat_wgs_two[dat_wgs_two$cisDist != '100', ]
dat_wgs_two = dat_wgs_two[dat_wgs_two$cisDist != '1000000', ]
dat_wgs_two = dat_wgs_two[dat_wgs_two$cisDist != '100000', ]
dat_wgs_two$cisDist = factor(dat_wgs_two$cisDist, levels = c("0", "300", "500", "1000", "10000", "100000"))

g = ggplot(data = dat_wgs_two) + 
  geom_bar(aes(x=cisDist, y=QTLs, fill=tool), stat = 'identity', position = "dodge") + 
  facet_wrap(~peak_calling) + 
  ylab('Number of QTLs at FDR < 0.05') + 
  xlab('Window size around Peak regions to include variants') + 
  ggtitle('Results from matrixEQTL using true genotype') + 
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
  scale_fill_brewer(palette = "Paired")

ggsave(paste0('results/matrixEQTL_and_fastQTL_results_true_genotype.png'), height = 2.5, width = 5)




dat_evaluation = read.table('data/Evaluation_metrics/fastQTL_results_variants_called_summary.txt', sep='\t', header = T, stringsAsFactors = F)
dat_evaluation$Recall = dat_evaluation$True_hits / dat_evaluation$True_QTLs
dat_evaluation$Precision = dat_evaluation$True_hits / dat_evaluation$QTLs
dat_evaluation$cisDist = factor(dat_evaluation$cisDist)
dat_evaluation$method = factor(dat_evaluation$method, levels = c("VCF_files", "Imputation", "Integration"),
                               labels = c("Genotype caller", "Imputation", "Integrated"))

g = ggplot(data = dat_evaluation) + 
  geom_bar(aes(x=cisDist, y=Recall, fill=method), stat = 'identity', position = "dodge") + 
  facet_wrap(~peak_calling) + 
  ylab('Recall') + 
  xlab('Window size around Peak middle point to include variants') + 
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
  scale_fill_viridis(discrete = TRUE, option = "D")

ggsave(paste0('results/fastQTL_results_Recall.png'), height = 2.5, width = 5)



g = ggplot(data = dat_evaluation) + 
  geom_bar(aes(x=cisDist, y=Precision, fill=method), stat = 'identity', position = "dodge") + 
  facet_wrap(~peak_calling) + 
  ylab('Precision') + 
  xlab('Window size around Peak middle point to include variants') + 
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
  scale_fill_viridis(discrete = TRUE, option = "D")

ggsave(paste0('results/fastQTL_results_Precision.png'), height = 2.5, width = 5)



files = list.files('data/Evaluation_metrics/')
files = files[grepl('topVariant.txt', files)]
dat = NULL
for(f in files){
  print(f)
  dt = read.table(paste0('data/Evaluation_metrics/', f), sep=',', header=T)
  dt$peak_calling = gsub('_minDP3_topVariant.txt', '', f)
  dat = rbind(dat, dt)
}

dat = dat[dat$peak_calling == 'MACS2', ]
dat$Window = factor(dat$Window, levels = c("300", "500", "1000", "10000", "100000", "1000000"))
#dat$peak_calling = factor(dat$peak_calling, levels = c("MACS2", "MACS2_combined", "Genrich", "Genrich_combined"),
#                                labels = c("MACS2 (call-merge)", "MACS2 (merge-call)", "Genrich (call-merge)", "Genrich (merge-call)"))
dat$method = factor(dat$method, levels = c("VCF_files", "Imputation", "Integration"),
                               labels = c("Genotype caller", "Imputation", "Integrated"))


g = ggplot(data = dat[(dat$THR == 1) | (dat$THR == 0.5)  , ]) + 
  geom_bar(aes(x=Window, y=Value, fill=method), stat = 'identity', position = "dodge") + 
  facet_wrap(~Metric + THR) + 
  xlab('Window size around Peak middle point to include variants') + 
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
  scale_fill_viridis(discrete = TRUE, option = "D")

ggsave(paste0('results/fastQTL_results_Precision.png'), height = 4, width = 6)
