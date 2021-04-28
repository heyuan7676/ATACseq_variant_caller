library(readr)
library(dplyr)

sample_structure = read.table('SraRunTables/correlation_Spearman.csv', sep='\t', header=T)
rownames(sample_structure) = sample_structure[, 1]
sample_structure = sample_structure[, rownames(sample_structure)]

samples = rownames(sample_structure)
#samples = samples[samples %in% colnames(vcf_dat)]
#samples = samples[samples %in% colnames(peak_dat)]
#sample_structure = sample_structure[samples, ]
#sample_structure = sample_structure[, samples]

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
		if( sum(sample_structure[s, samples_belong_to_this_donor] > 0.6) >= length(samples_belong_to_this_donor) / 2){
			sample_levels[s] = sample_levels[names(sample_levels)[potential_same_donor]] 
		}else{
			sample_group = sample_group + 1
			#print(sample_group)
			sample_levels[s] = paste0('S', sample_group)
		}
	}
}

samples = read.table('/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/merged_vcf_files/samples.txt', header = F, stringsAsFactors = F)
for(s in samples$V1){
	if(!(s %in% names(sample_levels))){
		sample_group = sample_group + 1
		sample_levels[s] = paste0('S', sample_group)
	}
}
write.table(data.frame(sample_levels), 'SraRunTables/correlation_Spearman_sample_structure.csv', sep='\t', quote =F)




