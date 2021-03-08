setwd('/Users/Yuan/Desktop/Variant_calling/')

library(ggplot2)
library(reshape2)

readin <-function(filename, method){
  dat = read.table(filename, sep=',', header=T, stringsAsFactors = F)
  dat$group = method
  return(dat)
}


plot_metric <- function(dat, savefile){
  dat$group = factor(dat$group,  levels = c("Genrich","MACS2", "Genrich_combined", "MACS2_combined"))
  dat$THR = as.factor(dat$THR)
  
  for(metric in c("Precision", "Recall")){
    df_plot = dat[dat$Metric == metric, ]
    g = ggplot(data = df_plot) + 
      geom_bar(aes(x=THR, y=Value, fill=method), stat = 'identity', position = "dodge") + 
      facet_wrap(~ Window + group, ncol = 4) + 
      ylab(metric) + 
      xlab('Criteria of same top QTL (R2 >= )') + 
      theme_bw()+
      theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
      scale_fill_viridis(discrete = TRUE, option = "D")
    
    ggsave(paste0('results/',savefile,'_', metric, '.png'), height = 5, width = 9)
  }

  df_plot$Window = as.factor(df_plot$Window)
  g = ggplot(data = df_plot) + 
    geom_bar(aes(x=Window, y=N_called, fill=method), stat = 'identity', position = "dodge") + 
    facet_wrap(~ group, ncol = 2) + 
    ylab("Number of QTLs at FDR < 0.05") + 
    xlab('Window size') + 
    theme_bw()+
    theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
    scale_fill_viridis(discrete = TRUE, option = "D")
  
  ggsave(paste0('results/',savefile,'_Number_QTLs.png'), height = 2.5, width = 5)
  
}

obtain_dat <- function(minDP = 2, datadir = 'data/Evaluation_metrics/'){
  
  dat1 = readin(paste0(datadir, 'Genrich_minDP', as.character(minDP),'_topVariant.txt'), 'Genrich')
  dat2 = readin(paste0(datadir, 'macs2_minDP', as.character(minDP),'_topVariant.txt'), 'MACS2')
  dat3 = readin(paste0(datadir, 'MACS2_combined_minDP', as.character(minDP),'_topVariant.txt'), 'MACS2_combined')
  dat4 = readin(paste0(datadir, 'Genrich_combined_minDP', as.character(minDP),'_topVariant.txt'), 'Genrich_combined')
  
  dat = rbind(dat1, dat2, dat3, dat4)
  #dat = dat[dat$THR == 0, ]
  dat$minDP = minDP
  
  return(dat)
}

dat = obtain_dat(minDP =3)

plot_metric(dat[dat$Standard == 'R2_threshold', ], paste0('QTL_metrics_THR0_Pvalue'))



dat = dat[dat$Standard == 'R2_threshold', ]
dat$group = factor(dat$group,  levels = c("Genrich","MACS2", "Genrich_combined", "MACS2_combined"))

dat_R05 = dat[dat$THR == 0.5, ]
for(metric in c("Precision", "Recall")){
  df_plot = dat_R05[dat_R05$Metric == metric, ]
  g = ggplot(data = df_plot) + 
    geom_bar(aes(x=group, y=Value, fill=method), stat = 'identity', position = "dodge") + 
    facet_wrap(~ Window) + 
    ylab(metric) + 
    xlab('') + 
    theme_bw()+
    theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
    scale_fill_viridis(discrete = TRUE, option = "D")
  ggsave(paste0('results/QTLs_R0.5_',metric,'_Number_QTLs.png'), height = 2.5, width = 5)
}

dat_R1 = dat[dat$THR == 1, ]
for(metric in c("Precision", "Recall")){
  df_plot = dat_R1[dat_R1$Metric == metric, ]
  g = ggplot(data = df_plot) + 
    geom_bar(aes(x=group, y=Value, fill=method), stat = 'identity', position = "dodge") + 
    facet_wrap(~ Window) + 
    ylab(metric) + 
    xlab('') + 
    theme_bw()+
    theme(axis.text.x  = element_text(angle = 45, hjust = 1))  + 
    scale_fill_viridis(discrete = TRUE, option = "D")
  ggsave(paste0('results/QTLs_R1_',metric,'_Number_QTLs.png'), height = 2.5, width = 5)
}

