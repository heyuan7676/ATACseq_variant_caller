setwd('/Users/Yuan/Desktop/Variant_calling/')

library(ggplot2)
library(reshape2)

readin <-function(filename, group_text){
  dat = read.table(filename, sep='\t', header=T, stringsAsFactors = F)
  dat$Precision = dat$Overlap / dat$Called_QTL
  dat$Recall = dat$Overlap / dat$True_QTL[1]
  
  dat$group = group_text
  
  return(dat)
}

plot_metric <- function(dat, col_to_use, text, savefile){
  df_plot = dat[,c(col_to_use, "GQ_threshold", "group")]
  df_plot$group = factor(df_plot$group, 
                         levels = c("Genrich_noWeights","Genrich_withWeights","MACS2_noWeights","MACS2_withWeights",
                                    "Genrich_combined_noWeights", "Genrich_combined_withWeights",
                                    "MACS2_combined_noWeights", "MACS2_combined_withWeights"))
  df_plot$GQ_threshold = factor(df_plot$GQ_threshold,
                                levels = c("PP>0", "PP>0.4", "PP>0.5", "PP=1"))
  
  g = ggplot(data = df_plot) + 
    geom_bar(aes_string(x='GQ_threshold', y=col_to_use, fill='group'), 
             stat = 'identity', position = "dodge") + 
    ylab(text) + 
    xlab('Variant genotype calling thresholds') + 
    theme_bw()+
    theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
    scale_fill_brewer(palette='Paired')
  
  ggsave(paste0('results/',savefile,'.png'), height = 3, width = 5)
}


plot_qtl_metrics <- function(window, datadir = 'data/Evaluation_metrics/'){
  
  dat1 = readin(paste0(datadir, 'Genrich_minDP2', window,'.txt'), 'Genrich_withWeights')
  dat2 = readin(paste0(datadir, 'Genrich_minDP2', window,'_noWeight.txt'), 'Genrich_noWeights')
  dat3 = readin(paste0(datadir, 'macs2_minDP2', window,'_noWeight.txt'), 'MACS2_noWeights')
  dat4 = readin(paste0(datadir, 'macs2_minDP2', window,'.txt'), 'MACS2_withWeights')
  dat5 = readin(paste0(datadir, 'macs2_combined_minDP2', window,'.txt'), 'MACS2_combined_withWeights')
  dat6 = readin(paste0(datadir, 'macs2_combined_minDP2', window,'_noWeight.txt'), 'MACS2_combined_noWeights')
  dat7 = readin(paste0(datadir, 'Genrich_combined_minDP2', window,'.txt'), 'Genrich_combined_withWeights')
  dat8 = readin(paste0(datadir, 'Genrich_combined_minDP2', window,'_noWeight.txt'), 'Genrich_combined_noWeights')
  
  
  dat = rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8)
  
  plot_metric(dat, 'Recall', 'Recall', paste0('Recall', window))
  plot_metric(dat,'Precision', 'Precision', paste0('Precision', window))
  plot_metric(dat, 'Called_QTL', 'Number of called QTL', paste0('Called_QTL', window))
}

plot_qtl_metrics('_0.0kb')
plot_qtl_metrics('_0.1kb')
plot_qtl_metrics('_1.0kb')


datadir = 'data/Evaluation_metrics/'
dat1 = readin(paste0(datadir, 'macs2_minDP2_0.0kb.txt'), '0.0kb')
dat2 = readin(paste0(datadir, 'macs2_minDP2_0.1kb.txt'), '0.1kb')
dat3 = readin(paste0(datadir, 'macs2_minDP2_1.0kb.txt'), '1.0kb')

dat = rbind(dat1, dat2, dat3)
dat = dat[dat$GQ_threshold == 'PP>0', ]


plot_metric_locs <- function(dat, col_to_use, text, savefile){
  df_plot = dat[,c(col_to_use, "GQ_threshold", "group")]
  df_plot$group = factor(df_plot$group, 
                         levels = c("0.0kb", '0.1kb', '1.0kb'))
  g = ggplot(data = df_plot) + 
    geom_bar(aes_string(x='group', y=col_to_use, fill='group'), 
             stat = 'identity', position = "dodge") + 
    ylab(text) + 
    xlab('Window size around the peaks to test for variants') + 
    theme_bw()+
    theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
    scale_fill_brewer(palette='Set1')
  
  ggsave(paste0('results/',savefile,'.png'), height = 3, width = 3)
}


plot_metric_locs(dat, 'Recall', 'Recall', paste0('Recall_different_locs'))
plot_metric_locs(dat, 'Precision', 'Precision', paste0('Precision_different_locs'))
plot_metric_locs(dat, 'Called_QTL', 'Number of called QTLs', paste0('Called_QTL_different_locs'))











plot_overlapping_QTLs_basedOn_pvalues <- function(peak_calling, minDP, overall_recall = F){
  dat = read.table(paste0('data/Evaluation_metrics/',peak_calling,'_minDP',minDP,'_withImputation_noWeight_topVariant.txt'), sep=',', header = T)

  #dat = dat[dat$THR == 0, ]
  
  if(overall_recall){
    dat$recall_called = dat$N_called * dat$recall_called / dat$N_true[dim(dat)[1]]
    dat$recall_imputed = dat$N_imputed * dat$recall_imputed / dat$N_true[dim(dat)[1]]
  }
  
  col = c("Window", "THR", "R2","precision_called","precision_imputed", "recall_called", "recall_imputed")
  metrics_QTLs = melt(dat[,col], id.vars  = c('Window', 'THR', "R2"))
  
  metrics_QTLs$metric = "METRIC"
  metrics_QTLs[grep("precision", metrics_QTLs$variable), "metric"] = "Precision"
  metrics_QTLs[grep('recall', metrics_QTLs$variable), "metric"] = "Recall"
  
  metrics_QTLs$group = "GROUP"
  metrics_QTLs[grep("called", metrics_QTLs$variable), "group"] = "Genotype caller"
  metrics_QTLs[grep("imputed", metrics_QTLs$variable), "group"] = "Imputation"
  
  metrics_QTLs$THR = factor(as.character(metrics_QTLs$THR), levels = rev(as.character(unique(metrics_QTLs$THR))))
  metrics_QTLs$R2 = factor(as.character(metrics_QTLs$R2), levels = as.character(unique(metrics_QTLs$R2)))
  metrics_QTLs$Window = factor(as.character(metrics_QTLs$Window), levels = c("0", "1000", "10000","100000"), 
                               labels = c("QLTs within Peaks", "QTLs <= 1kb around peaks", "QTLs <= 10kb around peaks","QTLs <= 100kb around peaks"))

  g = ggplot() + 
    geom_point(data = metrics_QTLs, aes(x = THR, y=value, color = group, group = group)) + 
    geom_line(data = metrics_QTLs, aes(x = THR, y=value, color = group, group = group)) + 
    facet_grid(metric~Window, scales = "free") + 
    theme_bw() + 
    ylab("") + 
    xlab("P-value Threshold to define overlapping") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color=guide_legend(title=""))
  
  ggsave(paste0('results/QTL_topVariant_overlapping_', peak_calling, '_minDP', minDP, '.png'), height = 4, width = 10)
  
}
plot_overlapping_QTLs_basedOn_pvalues('macs2', 2)

