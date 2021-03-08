setwd('/Users/Yuan/Desktop/Variant_calling/')

library(ggplot2)
library(reshape2)


plot_peak_percentile_by_QTLs <- function(peak_calling, peak_calling_text, window){
  peak_dat = read.table(paste0('data/peak_inPercentile_QTLs/Peak_percentile_QTLs_',peak_calling,'_',as.character(window),'.txt'), 
                        sep='\t', header=T, stringsAsFactors = F)
  
  breaks <- quantile(peak_dat$median_peak_pv, seq(0,10)/10)
  tags <- paste0("Percentile ", seq(1, 10) * 10, "%")
  
  group_tags <- cut(peak_dat$median_peak_pv, 
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=TRUE, 
                    labels=tags)
  peak_dat$peak_groups <- factor(group_tags,  levels = tags, ordered = TRUE)
  
  breaks <- quantile(peak_dat$number_samples, seq(0,5)/5)
  tags <- paste0("Number of samples (Percentile ", seq(1, 5) * 20, "%)")
  group_tags <- cut(peak_dat$number_samples, 
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=TRUE, 
                    labels=tags)
  peak_dat$sample_groups <- factor(group_tags,  levels = tags, ordered = TRUE)
  
  
  number_QTLs = peak_dat %>% group_by(peak_groups, sample_groups) %>% summarise(value = sum(bh < 0.05) / sum(bh > -1))
  
  g = ggplot(data = number_QTLs, aes(x=peak_groups, y=value, fill = sample_groups)) + 
    geom_bar(stat = 'identity') + 
    facet_wrap(~sample_groups) + 
    xlab(paste0('Median p-value percentile of peaks called by ', peak_calling)) +
    ylab('Fraction of tested peaks that yield ca-QTLs') + 
    ggtitle(paste0(peak_calling_text, ' variants within ',as.character(window),' bp')) + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle = 90))
  
  ggsave(paste0('results/peak_percentile_QTLs_',peak_calling,'_',as.character(window), '.png'), height = 5, width = 10)
  
  return(number_QTLs)
  
}


d1 = plot_peak_percentile_by_QTLs('MACS2', 'MACS2 (call-merge)', 1000)
d2 = plot_peak_percentile_by_QTLs('Genrich', 'Genrich (call-merge)', 1000)



plot_peak_percentile_by_QTLs_v1 <- function(peak_calling, peak_calling_text, window){
  peak_dat = read.table(paste0('data/peak_inPercentile_QTLs/Peak_percentile_QTLs_',peak_calling,'_',as.character(window),'.txt'), 
                        sep='\t', header=T, stringsAsFactors = F)
  
  breaks <- quantile(peak_dat$peak_pv, seq(0,10)/10)
  tags <- paste0("Percentile ", seq(1, 10) * 10, "%")
  
  group_tags <- cut(peak_dat$peak_pv, 
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=TRUE, 
                    labels=tags)
  peak_dat$peak_groups <- factor(group_tags,  levels = tags, ordered = TRUE)
  
  
  number_QTLs = peak_dat %>% group_by(peak_groups) %>% summarise(value = sum(bh < 0.05) / sum(bh > -1))
  
  g = ggplot(data = number_QTLs, aes(x=peak_groups, y=value)) + 
    geom_bar(stat = 'identity') + 
    xlab(paste0('Median p-value percentile of peaks called by ', peak_calling_text)) +
    ylab('Fraction of tested peaks that yield ca-QTLs') + 
    ggtitle(paste0(peak_calling_text, ' variants within ',as.character(window),' bp')) + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle = 90))
  
  ggsave(paste0('results/peak_percentile_QTLs_',peak_calling,'_',as.character(window), '.png'), height = 5, width = 6)
  
  return(number_QTLs)
  
}


plot_peak_percentile_by_QTLs <- function(peak_calling, window){
  print(c(peak_calling, window))
  peak_dat = read.table(paste0('data/peak_inPercentile_QTLs/Peak_percentile_QTLs_',peak_calling,'_',as.character(window),'.txt'), 
                        sep='\t', header=T, stringsAsFactors = F)
  
  breaks = c(-log10(0.05), -log10(1 / 10 ^ seq(2,22,2)), max(peak_dat$median_peak_qv))
  
  group_tags <- cut(peak_dat$median_peak_qv, 
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=TRUE)
  peak_dat$peak_groups <- factor(group_tags,  
                                 levels = levels(group_tags), 
                                 #labels = as.character(c(0.05, (1 / 10 ^ seq(2,22,2)))),
                                 ordered = TRUE)
  
  
  peak_dat_combined = peak_dat
  
  breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 71, 72)
  group_tags <- cut(peak_dat_combined$number_samples, 
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=TRUE)
  peak_dat_combined$sample_groups <- factor(group_tags,  
                                            levels = levels(group_tags), 
                                            ordered = TRUE)
  
  
  number_QTLs = peak_dat_combined %>% group_by(peak_groups, sample_groups) %>% summarise(value = sum(bh > -1) )
  g = ggplot(data = number_QTLs, aes(x=peak_groups, y=value, fill = sample_groups)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    xlab(paste0('Median -log10(q-value) threshold of peaks')) +
    #ylab('Fraction of tested peaks\n that yield ca-QTLs') + 
    ggtitle(paste0('# tested peaks')) + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle = 90))
  ggsave(paste0('results/peak_percentile_numberTests_',peak_calling,'_',as.character(window), '_qv.png'), height = 3, width = 5)
  
  number_QTLs = peak_dat_combined %>% group_by(peak_groups, sample_groups) %>% summarise(value = sum(bh < 0.05) )
  g = ggplot(data = number_QTLs, aes(x=peak_groups, y=value, fill = sample_groups)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    xlab(paste0('Median -log10(q-value) threshold of peaks')) +
    #ylab('Fraction of tested peaks\n that yield ca-QTLs') + 
    ggtitle(paste0('# peaks with QTLs')) + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle = 90))
  ggsave(paste0('results/peak_percentile_numberQTLs_',peak_calling,'_',as.character(window), '_qv.png'), height = 3, width = 5)
  
  number_QTLs = peak_dat_combined %>% group_by(peak_groups, sample_groups) %>% summarise(value = sum(bh < 0.05) / sum(bh > -1))
  g = ggplot(data = number_QTLs, aes(x=peak_groups, y=value, fill = sample_groups)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    xlab(paste0('Median -log10(q-value) threshold of peaks')) +
    ylab('Fraction of tested peaks\n that yield ca-QTLs') + 
    ggtitle(paste0('Fraction of peaks with QTLs')) + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle = 90))
  ggsave(paste0('results/peak_percentile_QTLs_',peak_calling,'_',as.character(window), '_qv.png'), height = 3, width = 5)
}




window = 300
peak_calling = 'MACS2'
for(window in c(300, 500, 1000, 10000, 100000, 1000000)){
  plot_peak_percentile_by_QTLs(peak_calling, window)  
}

peak_calling = 'MACS2'
files = list.files('data/peak_inPercentile_QTLs/')
files = files[grepl(peak_calling, files)]
files = files[grepl('subset', files)]

dat = NULL
for(f in files){
  print(f)
  dt = read.table(paste0('data/peak_inPercentile_QTLs/', f), sep='\t', header=T)
  dat = rbind(dat, dt)
}

dat$median_peak_qvalue = factor(dat$median_peak_qvalue, 
                                levels = c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001), 
                                labels = c("0.05", "0.01", "0.005", "0.001", "0.0005", "0.0001", "0.00005", "0.00001", "0.000005", "0.000001"))
dat$window = factor(dat$window)

g = ggplot(data = dat, aes(x = median_peak_qvalue, y = QTL_number, color=window, group=window)) + 
  geom_point()  + 
  geom_line() + 
  facet_wrap(~window) +
  xlab('Median peak qvalue threshold') + 
  ylab('QTL number at FDR < 0.05') + 
  scale_color_viridis(discrete = TRUE, option = "D") + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = 90))
ggsave(paste0('results/peak_subset_QTLs_',peak_calling, '_', as.character(window), '.png'), height = 4, width = 6)



# Pi1
peak_calling = 'MACS2'
dat = read.table(paste0('data/Evaluation_metrics/',peak_calling,'_minDP3_topVariant_cisDist_pi1_called_repInReal.txt'), sep='\t', header = T)
dat$discovery = "Discovery: Called"

dat2 = read.table(paste0('data/Evaluation_metrics/',peak_calling,'_minDP3_topVariant_cisDist_pi1_real_repInCalled.txt'), sep='\t', header = T)
dat2$discovery = "Discovery: Real"

dat = rbind(dat, dat2)
dat$window = factor(dat$window)
dat$method = factor(dat$method, levels = c("VCF_files", "Imputation", "Integration"),
                    labels = c("Genotype Caller", "Imputation", "Integration"))

g = ggplot(data = dat, aes(x = window, y = pi1, fill=method)) + 
  facet_wrap(~discovery)+ 
  geom_bar(stat = "identity", position = 'dodge') + 
  xlab('Window to test variants') + 
  ylab('pi1 ') + 
  scale_fill_viridis(discrete = TRUE, option = "D") + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = 90))

ggsave(paste0('results/peak_subset_QTLs_',peak_calling, '_pi1.png'), height = 2.3, width = 6)


