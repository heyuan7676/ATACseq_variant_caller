setwd('/Users/Yuan/Desktop/Variant_calling/')

library(ggplot2)

### Peaks
compare_peak_set_plot_overlapping <- function(filename, first_column, title, savefile, use_basepair = F){
  dat = read.table(filename, sep='\t', header = T)
  print(dat)
  df_plot = cbind(first_column, 
                  c("Unique", "Overlap", "Unique", "Overlap"))
  df_plot = as.data.frame(df_plot)
  df_plot$V2 = factor(df_plot$V2, levels = c("Unique", "Overlap"))
  
  if(use_basepair){
    df_plot$V3 = c(dat[1,4] - dat[1,6], dat[1,6], dat[1,5]-dat[1,6], dat[1,6])
    ylabel = 'Number of basepair in peak regions'
  }else{
    df_plot$V3 = c(dat[1,1] - dat[1,3], dat[1,3], dat[1,2]-dat[1,3], dat[1,3]) 
    ylabel = 'Number of peaks'
  }
  print(df_plot)
  
  g = ggplot() + 
    geom_bar(aes(x=first_column, y=V3, fill=V2), 
             data=df_plot, 
             stat='identity',
             position = 'stack') + 
    ylab(ylabel) + 
    xlab('') + 
    theme_bw() + 
    theme(axis.text.x = element_blank())  + 
    scale_fill_brewer(palette = "Set1") + 
    #scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
    guides(fill=guide_legend(title="")) + 
    ggtitle(title) + 
    theme_bw()
  ggsave(savefile, g, width = 3.3, height = 3.8)
}


first_col = c("s1", "s1", "s2", "s2")
s1 = 'SRR10972301'
s2 = 'SRR10972302'
compare_peak_set_plot_overlapping(paste0('data/GEO/evaluate_peaks/Examine_peaks_results_macs2_',s1,'_',s2,'.txt'), 
                                  first_col, 
                                  title = paste(s1, s2),
                                  savefile = paste0("results/Examine_peaks_overlap_",s1,"_",s2, ".png"), 
                                  use_basepair = F)

s1 = 'SRR10972301'
s2 = 'SRR10972302'
compare_peak_set_plot_overlapping(paste0('data/GEO/evaluate_peaks/Examine_peaks_results_macs2_',s1,'_',s2,'.txt'), 
                                  first_col, 
                                  title = paste(s1, s2),
                                  savefile = paste0("results/Examine_peaks_overlap_",s1,"_",s2, ".png"), 
                                  use_basepair = F)

s1 = 'SRR10972303'
s2 = 'SRR10972304'
compare_peak_set_plot_overlapping(paste0('data/GEO/evaluate_peaks/Examine_peaks_results_macs2_',s1,'_',s2,'.txt'), 
                                  first_col, 
                                  title = paste(s1, s2),
                                  savefile = paste0("results/Examine_peaks_overlap_",s1,"_",s2, ".png"), 
                                  use_basepair = F)

s1 = 'SRR10972305'
s2 = 'SRR10972306'
compare_peak_set_plot_overlapping(paste0('data/GEO/evaluate_peaks/Examine_peaks_results_macs2_',s1,'_',s2,'.txt'), 
                                  first_col, 
                                  title = paste(s1, s2),
                                  savefile = paste0("results/Examine_peaks_overlap_",s1,"_",s2, ".png"), 
                                  use_basepair = F)

s1 = 'SRR10972307'
s2 = 'SRR10972308'
compare_peak_set_plot_overlapping(paste0('data/GEO/evaluate_peaks/Examine_peaks_results_macs2_',s1,'_',s2,'.txt'), 
                                  first_col, 
                                  title = paste(s1, s2),
                                  savefile = paste0("results/Examine_peaks_overlap_",s1,"_",s2, ".png"), 
                                  use_basepair = F)

s1 = 'SRR10972309'
s2 = 'SRR10972310'
compare_peak_set_plot_overlapping(paste0('data/GEO/evaluate_peaks/Examine_peaks_results_macs2_',s1,'_',s2,'.txt'), 
                                  first_col, 
                                  title = paste(s1, s2),
                                  savefile = paste0("results/Examine_peaks_overlap_",s1,"_",s2, ".png"), 
                                  use_basepair = F)

### Peaks - overall
### Variants
files = list.files('data/GEO/evaluate_peaks/')
files = files[grep('Examine_peaks_results_macs2', files)]
dat = NULL
for(f in files){
  dati = read.table(paste0('data/GEO/evaluate_peaks/', f), sep='\t', header =T, stringsAsFactors = F)
  dat = rbind(dat, dati)
}

overlapping_peaks = c(dat$N_peaks_s1_s2/dat$N_peaks_s1, dat$N_peaks_s1_s2/dat$N_peaks_s2)
plot_overlapping_fraction = data.frame( "overlapping_fraction"= overlapping_peaks)

g = ggplot() + 
  geom_boxplot(data = plot_overlapping_fraction, aes(x = 1, y = overlapping_peaks))  + 
  ylab('Overlapping fraction') + 
  xlab('') + 
  theme_bw() + 
  guides(fill=guide_legend(title="")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave('results/GEO_overlapping_peaks.png', height = 3, width = 1.5)


### Variants
files = list.files('data/GEO/performance/overlapping/')
files = files[grep('overlapping', files)]
dat = NULL
for(f in files){
  dati = read.table(paste0('data/GEO/performance/overlapping/', f), sep=',', header =T, stringsAsFactors = F)
  dat = rbind(dat, dati)
}

overlapping_variants = c(dat$consistent/dat$only_df1, dat$consistent/dat$only_df2)
plot_overlapping_fraction = data.frame("imputation"=c(dat$include_imputation, dat$include_imputation), "overlapping_fraction"= overlapping_variants)
plot_overlapping_fraction$imputation = factor(as.character(plot_overlapping_fraction$imputation),
                                              levels = c("0", "1"),
                                              labels = c("Genotype caller", "Imputation"))

g = ggplot() + 
  geom_boxplot(data = plot_overlapping_fraction, aes(x = imputation, y = overlapping_variants,fill= imputation))  + 
  ylab('Overlapping fraction') + 
  xlab('') + 
  theme_bw() + 
  guides(fill=guide_legend(title="")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave('results/GEO_overlapping_genotype.png', height = 3, width = 3.5)




overlapping_variants = c(dat$consistent/dat$df1_df2)
plot_overlapping_fraction = data.frame("imputation"=c(dat$include_imputation), "overlapping_fraction"= overlapping_variants)
plot_overlapping_fraction$imputation = factor(as.character(plot_overlapping_fraction$imputation),
                                              levels = c("0", "1"),
                                              labels = c("Genotype caller", "Imputation"))

g = ggplot() + 
  geom_boxplot(data = plot_overlapping_fraction, aes(x = imputation, y = overlapping_variants,fill= imputation))  + 
  ylab('Overlapping fraction') + 
  xlab('') + 
  theme_bw() + 
  guides(fill=guide_legend(title="")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave('results/GEO_overlapping_genotype_among_bothcalled.png', height = 3, width = 3.5)



