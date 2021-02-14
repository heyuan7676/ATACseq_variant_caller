

args = commandArgs(trailingOnly=TRUE)
peak_calling = args[1]
cisDis = args[2]


outdir = paste0('/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/fastQTL/', peak_calling)
filename=paste0(outdir, '/window', as.character(cisDis),'.fastq.permutation.results.txt')


dat = read.table(filename, head=F, stringsAsFactors=F)
colnames(dat) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope","ppval", "bpval")
dat = dat[complete.cases(dat), ]
dat$bh = p.adjust(dat$bpval, method="fdr")
dat = dat[order(dat$bh), ]


outfile = paste0(outdir, '/window', as.character(cisDis),'.fastq.permutation.results.BH.txt')
write.table(dat[, c('pid', 'sid', 'dist', 'npval', 'ppval','bpval', 'bh')], outfile, quote=F, row.names=F, col.names=T)

print(paste('For', peak_calling, 'window size =', as.character(cisDis),': there are', sum(dat$bh < 0.05), 'significant QTLs with FDR < 0.05 among', dim(dat)[1], 'tested peaks'))

