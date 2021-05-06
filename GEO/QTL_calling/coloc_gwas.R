library(coloc)
library(readr)
library(data.table)
library(dplyr)

# Get all pairs for eQTL genes
args<-commandArgs(TRUE)
tissue = args[1]


eQTL_coloc_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/coloc_eQTLs/'
eQTL_coloc_dat = read_delim(paste0(eQTL_coloc_dir, tissue, '_eQTL_caQTL_coloc.txt'), delim=' ')
eQTL_coloc_dat = eQTL_coloc_dat %>% filter(H4 >= 0.8)

# e-Genes that are close to e-Peaks
nSample = read.table('/work-zfs/abattle4/heyuan/Variant_calling/GEO/QTL_calling/GTEx_sample_n.txt', sep='\t', header = F, stringsAsFactors = F)
nSample = nSample[nSample$V1 == tissue, "V2"]

eQTL_dir = '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/'
dat_eQTL = read_delim(paste0(eQTL_dir, tissue, ".allpairs.txt"), delim='\t')
dat_eQTL = dat_eQTL %>% filter(maf > 0.05)

maf_dat = dat_eQTL[,c('variant_id', 'maf')]
dat_eQTL = dat_eQTL %>% filter( gene_id %in% unique(eQTL_coloc_dat$gene))

include_vector = rep(0, dim(maf_dat)[1])
for(c in unique(eQTL_coloc_dat$X1)){
	print(c)
	include_vector = include_vector + grepl(paste0(c, '_'), maf_dat$variant_id)
}
maf_dat = maf_dat[which(include_vector > 0), ]
maf_dat = data.frame(maf_dat[!duplicated(maf_dat$variant_id), ])
rownames(maf_dat) = maf_dat$variant_id

maf_dat = as.data.table(maf_dat, keep.rownames=T)


# for each peak, perform coloc analysis
compute_coloc <- function(gwas_df, eQTL_df, peakid, N_data1 = 500){
  coloc_result = NULL
  for(geneid in unique(eQTL_df$gene_id)){
      print(geneid)
      eQTL_df_geneid = eQTL_df %>% filter(gene_id == geneid)
      N_data2 = unique(gwas_df$sample_size)[1]

      dataset1 = list()
      dataset1[['pvalues']] = eQTL_df_geneid$pval_nominal
      dataset1[['N']] = N_data1
      dataset1[['MAF']] = eQTL_df_geneid$maf
      dataset1[['beta']] = eQTL_df_geneid$slope
      dataset1[['varbeta']] = (eQTL_df_geneid$slope_se) ** 2
      dataset1[['type']] = 'quant'
      dataset1[['snp']] = eQTL_df_geneid$variant_id
      
      dataset2 = list()
      dataset2[['pvalues']] = gwas_df$pvalue
      dataset2[['N']] = N_data2
      dataset2[['MAF']] = as.data.frame(gwas_df$maf)$ma
      dataset2[['beta']] = gwas_df$effect_size
      dataset2[['varbeta']] = (gwas_df$standard_error) ** 2
      dataset2[['type']] = 'quant'
      dataset2[['snp']] = gwas_df$gtex_variant_id 

      if(length(intersect(dataset1[['snp']], dataset2[['snp']]) ) == 0){
	coloc_result_i = c(0, -1, -1, -1, -1, -1)
      }else{
      	coloc_result_i = as.vector(coloc.abf(dataset1, dataset2)$summary)
      }
      coloc_result = rbind(coloc_result, c(geneid, coloc_result_i))
  }
  return(coloc_result)
}

gwas_dir = '/work-zfs/abattle4/parsana/gtex_trans_v8/data/gwas_haky/'
outdir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/coloc_eQTLs/'


filenames = list.files(gwas_dir)
filenames = filenames[grepl('txt.gz', filenames)]
result = NULL
for(gwas_file in filenames){
	print(gwas_file)
	gwas_trait = gsub('.txt.gz', '', gwas_file)
        gwas_dat = read_delim(paste0(gwas_dir, gwas_file), delim = '\t')
        gwas_dat = gwas_dat %>% filter(chromosome %in% unique(eQTL_coloc_dat$X1))

        if(!is.na(gwas_dat$frequency[1])){
                gwas_dat$maf = as.vector(sapply(gwas_dat$frequency, function(x) min(c(x, 1-x))))
        }else{
		setkey(maf_dat, rn)
		gwas_dat$maf = maf_dat[gwas_dat$gtex_variant_id, "maf"]
		gwas_dat = gwas_dat[which(complete.cases(gwas_dat$maf)), ]
        }

	result_i = NULL
	result_i = tryCatch(
	{	compute_coloc(gwas_dat, dat_eQTL, peakid, N_data1 = nSample)
	},
	error = function(){
		print('no available data')
	}
	warning=function(){
		print('warning')
	}	
	)

	if(!is.null(result_i)){
                result_i = data.frame(result_i)
                result_i$trait = gwas_trait

                colnames(result_i) = c("geneid", "n_snps", "H0", "H1", "H2", "H3", "H4", "Trait")
                result_i$H0 = as.numeric(as.character(result_i$H0))
                result_i$H1 = as.numeric(as.character(result_i$H1))
                result_i$H2 = as.numeric(as.character(result_i$H2))
                result_i$H3 = as.numeric(as.character(result_i$H3))
                result_i$H4 = as.numeric(as.character(result_i$H4))

                write.table(result_i, paste0(outdir, tissue, '_GWAS_STUDY_',gwas_trait, '_eQTL_caQTL_gwasQTL_coloc.txt'), quote = F, row.names=F)

                result = rbind(result, result_i)
	}

}



result = merge(result, eQTL_coloc_dat, by.x= 'geneid', by.y= 'gene')
result = result[, c("trait", "geneid", "peak", "H4.x", "H4.y", "X1", "X2", "X3")]
colnames(result) = c('trait', 'geneid', 'peak', 'H4_GWAS_eQTL', 'H4_eQTL_caQTL', "peak_chr", 'peak_start', 'peak_end')
result = result[rev(order(result$H4_GWAS_eQTL)), ]
write.table(result, paste0(outdir, tissue, '_eQTL_caQTL_gwasQTL_coloc.txt'), quote = F, row.names=F)

