library(readr)
library(dplyr)

derive_donor <- function(sample_structure){
  samples = colnames(sample_structure)
  
  sample_levels = c()
  sample_group = 1
  sample_levels[samples[1]] = paste0('S', sample_group)
  
  for(s in samples[2: length(samples)]){
    potential_same_donor = which.max(sample_structure[s, names(sample_levels)])
    if(length(potential_same_donor) == 0){
      sample_group = sample_group + 1
      sample_levels[s] = paste0('S', sample_group)
    }else{
      samples_belong_to_this_donor = which(sample_levels == sample_levels[potential_same_donor])
      if( sum(sample_structure[s, samples_belong_to_this_donor] > 0.5) >= length(samples_belong_to_this_donor) / 2){
        sample_levels[s] = sample_levels[names(sample_levels)[potential_same_donor]]
      }else{
        sample_group = sample_group + 1
        print(sample_group)
        sample_levels[s] = paste0('S', sample_group)
      }
    }
  }
  
  for(s in samples){
    if(!(s %in% names(sample_levels))){
      sample_group = sample_group + 1
      sample_levels[s] = paste0('S', sample_group)
    }
  }
  return(sample_levels)
  
}


args = commandArgs(trailingOnly=TRUE)
study = args[1]
sample_structure = read.table(paste0('correlation/correlation_Pearson_minDP5_octopus_', study, '.csv'), sep='\t', header=T)  
  rownames(sample_structure) = sample_structure[, 1]
  sample_structure = sample_structure[, rownames(sample_structure)]
  
  #fastq_batch = sapply(rownames(sample_structure), function(x) strsplit(x, 'S')[[1]][3])
  #sample_structure = sample_structure[which(fastq_batch == dp_depth), which(fastq_batch == dp_depth)]
  sample_structure[is.na(sample_structure)] = 1
  assigned_donors = derive_donor(sample_structure)
  names(assigned_donors) = paste0("S", sapply(names(assigned_donors), function(x) strsplit(x, 'S')[[1]][2]))
  assigned_donors = data.frame(assigned_donors)
  assigned_donors$assigned_donors = paste0(study, '-', assigned_donors$assigned_donors)
  write.table(assigned_donors, paste0('correlation/samples_assigned_', study, '.txt'), sep='\t', quote = F)




