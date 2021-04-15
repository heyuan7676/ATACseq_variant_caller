

result = NULL

for(peak_calling in c('MACS2', 'Genrich', 'MACS2/combined', 'Genrich/combined')){
	for(cisDis in c('100','300', '500', '1000', '10000', '100000', '1000000')){
		print(c(peak_calling, cisDis))
		outdir = paste0('/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/fastQTL/', peak_calling)
		outfile = paste0(outdir, '/window', as.character(cisDis),'.fastq.permutation.results.BH.txt')
		dat = read.table(outfile, sep=' ', header =T)
		result = rbind(result, c(peak_calling, cisDis, dim(dat)[1], sum(dat$bh < 0.05)))
	}
}

result = as.data.frame(result)
colnames(result) = c('peak_calling', 'cisDist', 'Tested_peaks', 'QTLs')
print(result)
write.table(result, 'fastQTL_results_summary.txt', sep='\t', quote=F, row.names=F, col.names=T)
