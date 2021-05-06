library(readr)
library(data.table)
library(dplyr)

qtl_filename = 'window1000000.fastq.permutation.results.BH'

# Get ca-QTL results
fastqtl_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/fastQTL/'
caQTL_sig = read_delim(paste0(fastqtl_dir, qtl_filename, '.txt'), delim=' ', col_names = T)
caQTL_sig = caQTL_sig %>% filter(bh < 0.05)

# Add MAF for the variants
maf.dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/merged_vcf_files/'
maf.dat.all = NULL
for(chromosome in seq(1,22)){
        maf.file = paste0(maf.dir, 'chr', chromosome,'.maf005.biallelic.recode.maf.txt')
        maf.dat = read_delim(maf.file, delim=' ', col_names=F)
        colnames(maf.dat) = c('sid', 'maf')
        maf.dat = maf.dat %>% filter(sid %in% caQTL_sig$sid)
        maf.dat$maf = as.numeric(maf.dat$maf)
	maf.dat = maf.dat[complete.cases(maf.dat), ]
	maf.dat.all = rbind(maf.dat.all, maf.dat)
}
maf.dat.all$maf = as.vector(sapply(maf.dat.all$maf, function(x) min(c(x, 1-x))))
caQTL_sig_withMAF = merge(caQTL_sig, maf.dat.all, by.x = 'sid', by.y='sid')

caQTL_sig_withMAF$maf_bin = cut(caQTL_sig_withMAF$maf, breaks = seq(0,10)/10)
caQTL_sig_withMAF$distance_bin = cut(abs(caQTL_sig_withMAF$dist), breaks = seq(0,1000000, 10000))
caQTL_sig_withMAF = as_tibble(caQTL_sig_withMAF)
write.table(caQTL_sig_withMAF, paste0(fastqtl_dir, qtl_filename, '.withMAF.txt'), row.names = F, quote = F)


# Add MAF for the variants
dat_caQTL = read_delim(paste0(fastqtl_dir, 'random_pool.withMAF.txt'), delim=' ', col_names = F)
dat_caQTL$maf_bin = cut(dat_caQTL$X6, breaks = seq(0,10)/10)
dat_caQTL$distance_bin = cut(abs(dat_caQTL$X3), breaks = seq(0,1000000, 10000))
dat_caQTL = as_tibble(dat_caQTL)

random_df = NULL
for(mafi in levels(caQTL_sig_withMAF$maf_bin)){
	ref = caQTL_sig_withMAF %>% filter(maf_bin == mafi)
	match = dat_caQTL %>% filter(maf_bin == mafi)
	print(c(mafi, dim(ref)[1]))
	for(disi in levels(ref$distance_bin)){
		ref_i = ref %>% filter(distance_bin == disi)
		match_i = match %>% filter(distance_bin == disi)
		random_i = match_i[sample(dim(match_i)[1], dim(ref_i)[1] * 3), ]
		random_df = rbind(random_df, random_i)
	}
}


random_df = data.frame(random_df)
random_df = random_df[, colnames(random_df)[1:6]]
colnames(random_df) = c("sid", "pid", "dist", "p-value", "slope", "maf")
write.table(random_df, paste0(fastqtl_dir, qtl_filename, '.randomMatched.txt'), row.names = F, quote = F)

