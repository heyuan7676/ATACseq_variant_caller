library(readr)
library(dplyr)
library(lme4)

args = commandArgs(trailingOnly=TRUE)
peak_calling = 'MACS2'
method = 'Imputation'
cisDis = 1000
chromosome = args[1]

root_dir = "/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis//ATAC_seq/alignment_bowtie/"


peakDir = paste0(root_dir, "/Peaks_", peak_calling, "/")
peakFn  = paste0(peakDir, "/peak_by_sample_matrix_RPKM_corrected_chromosome", chromosome, "_forFastQTL.bed.gz")
peak_dat = read_delim(peakFn, delim='\t')

VCF_dir = paste0(root_dir, method, "/minDP3/")
VCFfn   = paste0(VCF_dir, "/dosage_by_sample_matrix_chr", chromosome, ".recode.forFastQTL.vcf.gz")
vcf_dat = read_delim(VCFfn, delim='\t', comment = '##')


fastQTL_dir = paste0(root_dir, '/fastQTL/', method, '/', peak_calling, '/')
fastQTL_fn = paste0(fastQTL_dir, 'chr',chromosome,'.window', cisDis,'.fastq.results')
fastQTL_dat = read_delim(fastQTL_fn, delim=' ', col_names=F)


sample_structure = read.table('SraRunTables/correlation.csv', sep='\t', header=T)
rownames(sample_structure) = sample_structure[, 1]
sample_structure = sample_structure[, rownames(sample_structure)]

samples = rownames(sample_structure)
sample_levels = c()
sample_group = 1
sample_levels[samples[1]] = paste0('S', sample_group)

for(s in samples[1: length(samples)]){
	same_donor = which(sample_structure[s, names(sample_levels)] > 0.6)
	if(length(same_donor) > 0){
		sample_levels[s] = sample_levels[names(sample_levels)[same_donor[1]]]
	}else{
		sample_group = sample_group + 1
		print(sample_group)
		sample_levels[s] = paste0('S', sample_group)
	}
}


r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

compute_LLM <- function(peakid, peak_dat, vcf_dat, fastQTL_dat){
	snps = fastQTL_dat %>% filter(X1 == peakid) %>% pull(X2)
	y = peak_dat %>% filter(PEAK == peakid) %>% select(all_of(samples))
	x = vcf_dat %>% filter(ID %in% snps) %>% select(all_of(samples))

	lm_dat = as.data.frame(t(rbind(as.data.frame(y), as.data.frame(x))))
	colnames(lm_dat) = c('peak', snps)
	lm_dat$sample_group = sample_levels

	result = NULL
	for(s in snps){
		mixed.lmer <- lmer(peak ~ get(s) + (1|sample_group), data = lm_dat)
		r2 = r2.corr.mer(mixed.lmer)
		result = rbind(result, c(peakid, s, r2))
	}
	return(result)	
}

result = do.call(rbind, sapply(unique(fastQTL_dat$X1), function(x) compute_LLM(x, peak_dat, vcf_dat, fastQTL_dat)))
result = as.data.frame(result)
colnames(result) = c("peak", "snp", "R2")
result$R2 = as.numeric(paste(result$R2))

result = result[rev(order(result$R2)), ]

write.table(result, paste0(root_dir, 'mixLL/', 'chr',chromosome,'.window', cisDis,'.mixLL.r2.results.txt'), row.names=F, quote = F)




