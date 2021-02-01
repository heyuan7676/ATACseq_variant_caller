library(ggplot2)
library(viridis)
library(reshape2)
library(cowplot)


obtain_metric <- function(di, sample, pattern, group){
  correlation_dat = NULL
 
  column_list = c("Dosage", "Imputed_dosage","Y_random_forest",  "Y_linear_regression", "Y_logistic_regression", "Y_ordinal_regression")
 
  cori = c(sample, pattern, "Spearman")
  for(col in column_list){
    cori = c(cori, cor(di[,"True_GT"], di[,col], method = 'spearman'))
  }
  correlation_dat = rbind(correlation_dat, cori)
  
  cori = c(sample, pattern, "MSE")
  for(col in column_list){
    cori = c(cori, sum((di[,col] - di[,"True_GT"]) ** 2) /dim(di)[1])
  }
  correlation_dat = rbind(correlation_dat, cori)
  correlation_dat = as.data.frame(correlation_dat)
  correlation_dat$group = group
  return(correlation_dat)
  
}


### variants
readin_variants_df <- function(pattern, suffix){
  files = list.files(paste0(training_dir, '/', pattern, '/output'))
  files = files[grepl(suffix, files)]
  metrics_all = NULL
  metrics_by_GT = NULL
  metrics_by_dis = NULL
  
  for(filename in files){
    print(filename)
    dati = read.table(paste0(training_dir, '/', pattern, '/output/', filename), sep='\t', header = T, stringsAsFactors = F)
    sample = strsplit(filename, '_')[[1]][1]
    dati$sample = sample
    
    metrics_all = rbind(metrics_all, obtain_metric(dati, sample, pattern, "all"))
    
    for(gt in c(0,1,2)){
      metrics_by_GT = rbind(metrics_by_GT, obtain_metric(dati[dati$True_GT == gt, ], sample, pattern, as.character(gt))) 
    }
  
    for(distance in c(1e3, 1e4, 1e5)){
      metrics_by_dis = rbind(metrics_by_dis, obtain_metric(dati[dati$Distance_to_peaks <= log10(distance), ], sample, pattern, as.character(distance))) 
    }
    
    dati$True_GT = factor(dati$True_GT, levels = c(0, 1, 2), labels = c("AA", "AB", "BB"))
    dati = dati[sample(seq(1,dim(dati)[1]), 10000), ]
    g1 = ggplot(data = dati) + geom_boxplot(aes(x = True_GT, y = Y_random_forest)) + ylim(0, 2)
    g2 = ggplot(data = dati) + geom_boxplot(aes(x = True_GT, y = Y_linear_regression)) + ylim(0, 2)
    g3 = ggplot(data = dati) + geom_boxplot(aes(x = True_GT, y = Y_logistic_regression)) + ylim(0, 2)
    g4 = ggplot(data = dati) + geom_boxplot(aes(x = True_GT, y = Dosage)) + ylim(0, 2)
    g5 = ggplot(data = dati) + geom_boxplot(aes(x = True_GT, y = Imputed_dosage)) + ylim(0, 2)
    
    g = plot_grid(g4, g5, g1, g2, g3, nrow = 1)
    save_plot(paste0(training_dir, '/', pattern, '/output/', gsub(".txt", ".png", filename)), g, base_width = 10)
    
  }
  
  metrics_all = as.data.frame(metrics_all)
  metrics_by_GT = as.data.frame(metrics_by_GT)
  metrics_by_dis = as.data.frame(metrics_by_dis)
 
  colnames(metrics_all) = c("sample", "pattern", "Metric","GC", "IP", "RF", "LiR", "LgR", "OR", "Group")
  colnames(metrics_by_GT) = c("sample", "pattern", "Metric","GC", "IP", "RF", "LiR", "LgR", "OR", "Group")
  colnames(metrics_by_dis) = c("sample", "pattern", "Metric","GC", "IP", "RF",  "LiR", "LgR", "OR", "Group")
  
  return(list(metrics_all, metrics_by_GT, metrics_by_dis))
}



options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

minDP = as.numeric(args[1]) # 2
SUFFIX = args[2] # '_THR0.1.txt'

training_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model'
print(paste0('minDP', minDP))

saved_file = paste0('Performance/minDP', minDP, SUFFIX, '.RData')
if( file.exists(saved_file)){
  load(saved_file)
}else{
  resulti = readin_variants_df(paste0('minDP', minDP), SUFFIX)
  save(resulti, file = paste0('Performance/minDP', minDP, SUFFIX, '.RData'))
  print(resulti[[1]])
}


