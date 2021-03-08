setwd('/Users/Yuan/Desktop/Variant_calling/')

library(ggplot2)
library(reshape2)
library(VennDiagram)

dat = read.table('data/evaluate_peaks/Examine_peaks_results_genrich_macs2.txt', sep='\t', header =T)

# Peak overlap
overlap_peaks = dat[,c("N_peaks_macs2", "N_peaks_Genrich", "N_peaks_macs2_in_Genrich", "N_peaks_Genrich_in_macs2")]
colnames(overlap_peaks) = c("MACS2", "Genrich", "Overlap1", "Overlap2")

overlap_peaks$MACS2_only = overlap_peaks$MACS2 - overlap_peaks$Overlap1
overlap_peaks$Genrich_only = overlap_peaks$Genrich -overlap_peaks$Overlap2

df1 = overlap_peaks[,c("MACS2_only", "Overlap1")]
df1$sample = rownames(df1)
df1$group = "MACS2"
colnames(df1) = c("Unique", "Overlap", "Sample","Tool")
df1 = melt(df1, id.vars = c("Sample", "Tool"))

df2 = overlap_peaks[,c("Genrich_only", "Overlap2")]
df2$sample = rownames(df2)
df2$group = "Genrich"
colnames(df2) = c("Unique", "Overlap", "Sample","Tool")
df2 = melt(df2, id.vars = c('Sample', 'Tool'))

overlap_peaks_plot = rbind(df1, df2)


g = ggplot() + 
  geom_bar(aes(x=Sample, y=value, fill=variable), 
           data=overlap_peaks_plot, 
           stat='identity',
           position = 'stack') + 
  facet_wrap(~ Tool,nrow=2) + 
  ylab('Number of peaks') + 
  xlab('Samples') + 
  theme_bw() + 
  theme(axis.text.x = element_blank())  + 
  scale_fill_brewer(palette = "Set1")

ggsave("results/Examine_peaks_overlap_MACS2_Genrich.png", g)



# basepair overlap
overlap_peaks = dat[,c("bp_macs2", "bp_Genrich", "bp_overlap")]
colnames(overlap_peaks) = c("MACS2", "Genrich", "Overlap")
overlap_peaks$MACS2_only = overlap_peaks$MACS2 - overlap_peaks$Overlap
overlap_peaks$Genrich_only = overlap_peaks$Genrich -overlap_peaks$Overlap

df = overlap_peaks[,c("MACS2_only", "Genrich_only", "Overlap")]
df$sample = rownames(df)
df = melt(df, id.vars = "sample")

g = ggplot() + 
  geom_bar(aes(x=sample, y=value, fill=variable), 
           data=df, 
           stat='identity',
           position = 'stack') +
  ylab('Number of basepairs covered in peaks') + 
  xlab('Samples') + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_brewer(palette = "Blues")

ggsave("results/Examine_peaks_overlap_MACS2_Genrich_bp.png", g, height = 3)






####################################################################################################



compare_peak_set_plot <- function(filename, first_column, title, savefile, use_basepair = F){
  dat = read.table(filename, sep='\t', header = T)
  print(dat)
  df_plot = cbind(first_column, 
                  c("Unique", "Overlap", "Unique", "Overlap"))
  df_plot = as.data.frame(df_plot)
  df_plot$V2 = factor(df_plot$V2, levels = c("Unique", "Overlap"))
  
  if(use_basepair){
    df_plot$V3 = c(dat[1,5] - dat[1,7], dat[1,7], dat[1,6]-dat[1,7], dat[1,7])
    ylabel = 'Number of basepair in peak regions'
  }else{
    df_plot$V3 = c(dat[1,1] - dat[1,4], dat[1,4], dat[1,2]-dat[1,3], dat[1,3]) 
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




########## macs2 vs. genrich call-merge

file1 = 'data/evaluate_peaks/Examine_peaks_results_genrich_macs2_union.txt'
first_col = c("MACS2", "MACS2", "Genrich", "Genrich")

compare_peak_set_plot(file1, 
                         first_col, 
                         title = 'Call-merge',
                         savefile = "results/Examine_peaks_overlap_MACS2_Genrich_union.png")

compare_peak_set_plot(file1, 
                         first_col, 
                         title = 'Call-merge',
                         savefile = "results/Examine_peaks_overlap_MACS2_Genrich_union_bp.png", 
                         use_basepair = T)


########## macs2, call-merge vs. merge-call

file2 = 'data/evaluate_peaks/Examine_peaks_results_combined_macs2.txt'
first_col = c("Call_merge", "Call_merge", "Merge_call", "Merge_call")

compare_peak_set_plot(file2, 
                         first_col, 
                         title = 'MACS2',
                         savefile = "results/Examine_peaks_overlap_MACS2_combined.png")

compare_peak_set_plot(file2, 
                         first_col, 
                         title = 'MACS2',
                         savefile = "results/Examine_peaks_overlap_MACS2_combined_bp.png", 
                         use_basepair = T)



########## Genrich, call-merge vs. merge-call

file3 = 'data/evaluate_peaks/Examine_peaks_results_combined_genrich.txt'
first_col = c("Call_merge", "Call_merge", "Merge_call", "Merge_call")

compare_peak_set_plot(file3, 
                      first_col, 
                      title = 'Genrich',
                      savefile = "results/Examine_peaks_overlap_Genrich_combined.png")

compare_peak_set_plot(file3, 
                      first_col, 
                      title = 'Genrich',
                      savefile = "results/Examine_peaks_overlap_Genrich_combined_bp.png", 
                      use_basepair = T)



########## macs2 vs. Genrich merge-call

file4 = 'data/evaluate_peaks/Examine_peaks_results_genrich_macs2_combined.txt'
first_col = c("MACS2", "MACS2", "Genrich", "Genrich")

compare_peak_set_plot(file4, 
                      first_col, 
                      title = 'Merge-call',
                      savefile = "results/Examine_peaks_overlap_MACS2_Genrich_combined.png")

compare_peak_set_plot(file4, 
                      first_col, 
                      title = 'Merge-call',
                      savefile = "results/Examine_peaks_overlap_MACS2_Genrich_combined_bp.png", 
                      use_basepair = T)



