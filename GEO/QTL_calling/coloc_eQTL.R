library(coloc)
library(readr)
library(data.table)
library(dplyr)

# Get all pairs for eQTL genes
args<-commandArgs(TRUE)
tissue = args[1]

# restrict to SNPs that overlap with gwashits
gwas_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/coloc_gwasQTLs/'
gwas_hits_files = paste0(gwas_dir, 'global_', tissue, '_hits.txt')
gwas_hits = read.table(gwas_hits_files, header = F, stringsAsFactors = F)

# e-Genes that are close to e-Peaks
gene_close_to_peaks = read_delim('/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/Peaks_MACS2/genes_near_peaks.bed', delim='\t', col_names = F)

if(length(intersect(gwas_hits$V1, gene_close_to_peaks$X8)) == 0){
        print('No available eGenes')
        exit()
}

gene_close_to_peaks = gene_close_to_peaks %>% filter(X8 %in% gwas_hits$V1)


nSample = read.table('GTEx_sample_n.txt', sep='\t', header = F, stringsAsFactors = F)
nSample = nSample[nSample$V1 == tissue, "V2"]
eQTL_sig_dir = '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/'
dat_eGenes = fread(paste0(eQTL_sig_dir, tissue, ".v8.signif_variant_gene_pairs.txt"),select=c(2))
eGenes = unique(dat_eGenes$gene_id)
N0 = length(eGenes)

eQTL_dir = '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/'
dat_eQTL = read_delim(paste0(eQTL_dir, tissue, ".allpairs.txt"), delim='\t')
dat_eQTL = dat_eQTL[complete.cases(dat_eQTL), ]
dat_eQTL = dat_eQTL %>% filter(gene_id %in% gene_close_to_peaks$X8)
#dat_eQTL = dat_eQTL %>% filter(gene_id %in% eGenes)
#dat_eQTL = dat_eQTL %>% filter(abs(tss_distance) < 1000000)

# restrict to SNPs that overlap with gwashits

# Get ca-QTL results
fastqtl_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/fastQTL/'
caPeaks = read_delim(paste0(fastqtl_dir, 'window1000000.fastq.permutation.results.BH.txt'), delim=' ', col_names = T)
caPeaks = unique(caPeaks %>% filter(bh < 0.05) %>% pull(pid))
N1 =length(caPeaks)
dat_caQTL = read_delim(paste0(fastqtl_dir, 'window1000000.fastq.all.results.txt'), delim=' ', col_names = F)
dat_caQTL = dat_caQTL %>% filter(X1 %in% caPeaks)


gene_close_to_peaks = gene_close_to_peaks %>% filter(X8 %in% dat_eQTL$gene_id) %>% filter(X4 %in% dat_caQTL$X1)
gene_pool = unique(gene_close_to_peaks$X8)
dat_eQTL = dat_eQTL %>% filter(gene_id %in% gene_pool)
if(dim(dat_eQTL) == 0){
	print('No available eGenes')
	exit()
}
print(paste0(length(gene_pool), "/", N0," e-Genes located within 100kb of ca-QTL peaks"))

peak_pool = unique(gene_close_to_peaks$X4)
dat_caQTL = dat_caQTL %>% filter(X1 %in% peak_pool)
print(paste0(length(peak_pool),"/", N1, " e-Peaks located within 100kb of e-Genes"))


# Add MAF for the variants
maf.dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/merged_vcf_files/'
maf.dat.all = NULL
for(chromosome in seq(1,22)){
        maf.file = paste0(maf.dir, 'chr', chromosome,'.maf005.biallelic.recode.maf.txt')
        maf.dat = read_delim(maf.file, delim='  ', col_names=F)
        colnames(maf.dat) = c('sid', 'maf')
        maf.dat = maf.dat %>% filter(sid %in% dat_caQTL$X2)
        maf.dat$maf = as.numeric(maf.dat$maf)
	maf.dat = maf.dat[complete.cases(maf.dat), ]
	maf.dat.all = rbind(maf.dat.all, maf.dat)
}
maf.dat.all$maf = as.vector(sapply(maf.dat.all$maf, function(x) min(c(x, 1-x))))


#dat_caQTL_withMAF = merge(dat_caQTL, maf.dat.all, by.x = 'X2', by.y='sid')
#pv = apply(dat_caQTL_withMAF, 1, function(x) ifelse(x[4]>0, x[3], 1-x[3]))
#z_score = qnorm(pv)
#dat_caQTL_withMAF$se = dat_caQTL_withMAF$X5 / z_score


# for each peak, perform coloc analysis
compute_coloc <- function(caQTL_df, eQTL_df, peakid, N_data1 = 500, N_data2 = 838){
  caQTL_dat_for_coloc = caQTL_df %>% filter(X1 == peakid)
 
  genes_to_test = gene_close_to_peaks %>% filter(X4 == peakid) %>% pull(X8)
  eQTL_dat_for_coloc = eQTL_df %>% filter(gene_id %in% genes_to_test)
  geneids = unique(eQTL_dat_for_coloc %>% pull(gene_id))
  
  coloc_result = NULL
  
  for(geneid in geneids){
    eQTL_dat_for_coloc_genei = eQTL_dat_for_coloc %>% filter(gene_id == geneid)
    if(dim(eQTL_dat_for_coloc_genei)[1] > 0){
      dataset1 = list()
      dataset1[['pvalues']] = eQTL_dat_for_coloc_genei$pval_nominal
      dataset1[['N']] = N_data1
      dataset1[['MAF']] = eQTL_dat_for_coloc_genei$maf
      dataset1[['beta']] = eQTL_dat_for_coloc_genei$slope
      dataset1[['varbeta']] = (eQTL_dat_for_coloc_genei$slope_se) ** 2
      dataset1[['type']] = 'quant'
      dataset1[['snp']] = eQTL_dat_for_coloc_genei$variant_id
      
      dataset2 = list()
      dataset2[['pvalues']] = caQTL_dat_for_coloc$X4
      dataset2[['N']] = N_data2
      dataset2[['MAF']] = caQTL_dat_for_coloc$maf
      dataset2[['beta']] = caQTL_dat_for_coloc$X5
      #dataset2[['varbeta']] = (caQTL_dat_for_coloc$se) ** 2
      dataset2[['type']] = 'quant'
      dataset2[['snp']] = caQTL_dat_for_coloc$X2 

      print(c(geneid, peakid, min(eQTL_dat_for_coloc_genei$pval_nominal), min(caQTL_dat_for_coloc$X4)))

      if(length(intersect(dataset1[['snp']], dataset2[['snp']]) ) == 0){
	coloc_result_i = c(0, -1, -1, -1, -1, -1)
      }else{
      	coloc_result_i = as.vector(coloc.abf(dataset1, dataset2)$summary)
      }
      coloc_result_i = c(peakid, geneid, coloc_result_i)
      coloc_result = rbind(coloc_result, coloc_result_i)
      }
  }
  return(coloc_result)
}


dat_caQTL_withMAF$X2 = paste0('chr', dat_caQTL_withMAF$X2)
dat_eQTL$variant_id = as.vector(sapply(dat_eQTL$variant_id , function(x) paste(strsplit(x, '_')[[1]][1], strsplit(x, '_')[[1]][2], sep='_')))
peak_pool = intersect(peak_pool, dat_caQTL_withMAF$X1)
result = NULL
for(peakid in peak_pool){
  print(peakid)
  resi = compute_coloc(dat_caQTL_withMAF, dat_eQTL, peakid, N_data1 = nSample)
  if(length(resi) > 0){
      result = rbind(resi, result)
  }
}

result = data.frame(result)
colnames(result) = c("peak", "gene", "n_snps", "H0", "H1", "H2", "H3", "H4")
result$H0 = as.numeric(as.character(result$H0))
result$H1 = as.numeric(as.character(result$H1))
result$H2 = as.numeric(as.character(result$H2))
result$H3 = as.numeric(as.character(result$H3))
result$H4 = as.numeric(as.character(result$H4)) 

peak_loc = read_delim('/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/Peaks_MACS2/union-peaks.bed', delim = '\t', col_names = F)
result = merge(result, peak_loc, by.x = 'peak', by.y = 'X4')


result = result[rev(order(result$H4)), ]
outdir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/coloc_eQTLs/'
write.table(result, paste0(outdir, tissue, '_eQTL_caQTL_coloc.txt'), quote = F, row.names=F)

