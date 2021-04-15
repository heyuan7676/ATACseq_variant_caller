

result = NULL

for(peak_calling in c('MACS2', 'Genrich', 'MACS2/combined', 'Genrich/combined')){
	for(cisDis in c('500', '1000', '10000', '100000', '1000000')){
		print(c(peak_calling, cisDis))

		# true hits
		outdir = paste0('/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/fastQTL/', peak_calling)
		outfile = paste0(outdir, '/window', as.character(cisDis),'.fastq.permutation.results.BH.txt')
		dat = read.table(outfile, sep=' ', header =T)

		# called hits
		for(method in c('VCF_files', 'Imputation', 'Integration')){
			outdir = paste0('/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/fastQTL/', method, '/', peak_calling)
                	outfile = paste0(outdir, '/window', as.character(cisDis),'.fastq.permutation.results.BH.txt')
                	dat_called = read.table(outfile, sep=' ', header =T)
			d1 = dat_called[dat_called$bh < 0.05, ]
			d2 = dat[dat$bh < 0.05, ]
			result = rbind(result, c(peak_calling, cisDis, method, dim(dat_called)[1], sum(dat_called$bh < 0.05),  dim(merge(d1, d2, by = c('pid', 'sid')))[1], sum(dat$bh < 0.05)))
		}
	}
}

result = as.data.frame(result)
colnames(result) = c('peak_calling', 'cisDist', 'method', 'Tested_peaks', 'QTLs', 'True_hits', 'True_QTLs')
print(result)
write.table(result, 'fastQTL_results_variants_called_summary.txt', sep='\t', quote=F, row.names=F, col.names=T)
