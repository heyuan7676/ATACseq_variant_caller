# Matrix eQTL by Andrey A. Shabalin

# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)


args <- commandArgs(trailingOnly = TRUE)
CHROMOSOME = args[1]
method = args[2]
cisDist = as.numeric(args[3])

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/'

# Genotype file name
snp.dir = paste(root_dir, 'Dosage_for_QTL/minDP3/', method,'/', sep="") 
SNP_file_name = paste(snp.dir, "called_genotypes_chromosome", CHROMOSOME,".txt", sep="");
snps_location_file_name = paste(snp.dir, "called_genotypes_chromosome", CHROMOSOME,"_loc.bed", sep="");


# Gene expression file name
peak.dir = paste(root_dir, 'Peaks_MACS2/', sep="")
expression_file_name = paste(peak.dir, "peak_by_sample_matrix_RPKM_corrected_chromosome", CHROMOSOME,"_values.txt", sep="");
gene_location_file_name = paste(peak.dir, "peak_by_sample_matrix_RPKM_corrected_chromosome", CHROMOSOME,"_loc.bed", sep="");

# Covariates file name
# Set to character() for no covariates

# Output file name
qtl.dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/QTLs_MACS2/'
output_file_name_cis = paste(qtl.dir, "matrixeQTL_chromosome", CHROMOSOME,"_minDP3_cisDist_", as.character(cisDist/1000),"kb_", method,".txt", sep="");
output_file_name_tra = paste(qtl.dir, "matrixeQTL_chromosome", CHROMOSOME,"_minDP3_cisDist_", as.character(cisDist/1000),"kb_trans_", method, ".txt", sep="");

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

idx = sapply(colnames(snps), function(x) grep(x, colnames(gene)))
gene$ColumnSubsample(idx);

## Load covariates

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps,
gene = gene,
#cvrt = cvrt,
useModel = useModel,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
sum(me$cis$eqtls$FDR < 0.05)

