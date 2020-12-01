import sys
import time
import os
import pdb
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from statsmodels.regression import linear_model as sm

sys.path.append('/work-zfs/abattle4/heyuan/Variant_calling/scripts/QTL')
from prepare_data_matrix import *
from QTL_calling import *

def read_in_peaks(PEAK_dir, chromosome):
    peak_dat = pd.read_csv('%s/peak_by_sample_matrix_RPKM_corrected_chromosome%d.txt' % (PEAK_dir, chromosome), sep='\t')
    samples = [x for x in peak_dat.columns if x.startswith('HG')]
    return [peak_dat, samples]


if __name__ == "__main__":
    ROOT_DIR = sys.argv[1]  # /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/
    CHROMOSOME = int(sys.argv[2])
    WINDOW = int(sys.argv[3])
    peak_calling = sys.argv[4]
    imputation = True

    alignment_dir = 'alignment_bowtie' # alignment_subsample_0.5
    GT_subDir = 'minDP2'

    root_dir = '%s/%s' % (ROOT_DIR, alignment_dir)
    VCF_dir = '%s/VCF_files' % root_dir
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/maf005'

    if peak_calling == 'macs2':
    	PEAK_dir = '%s/Peaks' % root_dir
        QTL_dir = '%s/QTLs/%s' % (root_dir, GT_subDir)

    if peak_calling == 'macs2_combined':
        PEAK_dir = '%s/Peaks/combined' % root_dir
        QTL_dir = '%s/QTLs_combined/%s' % (root_dir, GT_subDir)

    elif peak_calling == 'Genrich':
        PEAK_dir = '%s/Peaks_Genrich' % root_dir
        QTL_dir = '%s/QTLs_Genrich/%s' % (root_dir, GT_subDir)

    elif peak_calling == 'Genrich_combined':
        PEAK_dir = '%s/Peaks_Genrich/combined' % root_dir
        QTL_dir = '%s/QTLs_Genrich_combined/%s' % (root_dir, GT_subDir)


    if imputation:
        QTL_dir = '%s/Imputation' % QTL_dir

    if not os.path.exists(QTL_dir):
        os.makedirs(QTL_dir)


    print('Call ca-QTLs for chromosome %d with window = %skb' % (CHROMOSOME, str(WINDOW/1000.0)))

    ## read in data
    [PEAK_DAT, SAMPLES_Peaks] = read_in_peaks(PEAK_dir, CHROMOSOME)

    ## align the samples with samples from WGS
    WGS_fn = '%s/%s.genotypes.tsv' % (WGS_dir, '1k_genome_chr22')
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t', nrows = 10)
    SAMPLES_WGS  = [x for x in WGS_result.columns if x.startswith('HG')]
    SAMPLES = list(np.intersect1d(SAMPLES_Peaks, SAMPLES_WGS))
    print('%d samples are used for QTL analysis\n' % len(SAMPLES))

    print('Read in genotype data from WGS...')
    [WGS_DAT, _] = read_in_WGS_GT('1k_genome_chr%d' % CHROMOSOME, SAMPLES, WGS_dir)
    print('done\n')

    ## read in data
    print('Read in genotype data from ATAC-seq reads...')
    Genotype_dir = '%s/Called_GT/%s' % (root_dir, GT_subDir)
    GT_DAT = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, snps = np.array(WGS_DAT['CHR_POS']))
    print('done\n')

    if imputation:
        print('Read in genotype data from imputation...')
        Genotype_dir = '%s/Imputation' % root_dir
        GT_DAT_Imputed = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, snps = np.array(WGS_DAT['CHR_POS']))
  	print('done\n')

	only_atac_reads = set(GT_DAT['CHR_POS']) - set(GT_DAT_Imputed['CHR_POS'])
	only_imputed = set(GT_DAT_Imputed['CHR_POS']) - set(GT_DAT['CHR_POS'])
	both = np.intersect1d(GT_DAT['CHR_POS'], GT_DAT_Imputed['CHR_POS'])

	GT_DAT_atac = GT_DAT.set_index('CHR_POS').loc[only_atac_reads]
	GT_DAT_imputed = GT_DAT_Imputed.set_index('CHR_POS').loc[only_imputed]

	# deal with in-consistent calls from the two sources
	GT_DAT_both = GT_DAT.set_index('CHR_POS').loc[both]
	imputed_both = GT_DAT_Imputed.set_index('CHR_POS').loc[both]
	genotype_both = np.array(imputed_both[SAMPLES])

	# remove genotypes that: 1). called by no method a heterozygous site; 2). called differently by imputed data and the ATAC-seq reads
	discard_snps = np.where((np.array(GT_DAT_both[SAMPLES]) != np.array(imputed_both[SAMPLES])) * (np.array(GT_DAT_both[SAMPLES]) != 1) * (np.array(imputed_both[SAMPLES]) != 1))
	genotype_both[discard_snps] = -1
	
	# use heterozygous calls from ATAC-seq reads
	ht_calls_from_reads = np.where((np.array(imputed_both[SAMPLES]) != 1) & (np.array(GT_DAT_both[SAMPLES]) == 1))
	genotype_both[ht_calls_from_reads] = 1 

	genotype_both = pd.DataFrame(genotype_both)
	genotype_both.columns = SAMPLES
	genotype_both.index = GT_DAT_both.index
	genotype_both['CHR'] = GT_DAT_both['CHR']
	genotype_both['POS'] = GT_DAT_both['POS']
	genotype_both = genotype_both[GT_DAT_imputed.columns]
	assert np.sum(np.array(GT_DAT_both)[np.where(np.array(genotype_both)==-1)] == 1)  == 0
	assert np.sum(np.array(imputed_both)[np.where(np.array(genotype_both)==-1)] == 1)  == 0
	
	GT_MERGED_DAT = GT_DAT_atac.append(GT_DAT_imputed).append(genotype_both)
	GT_MERGED_DAT['CHR_POS'] = GT_MERGED_DAT.index

	print("Obtained %d variants from ATAC-seq reads only, %d variants from imputation only, and %d variants from both sources" % (len(GT_DAT_atac), len(GT_DAT_imputed), len(GT_DAT_both)))	
	print("Obtain %d variants in total; discard %d(%.2f) snps' genotype that have different genotype information from atac-seq reads and imputation, and that are not heterozygous" % (len(GT_MERGED_DAT), len(discard_snps[0]), float(len(discard_snps[0])) / len(GT_MERGED_DAT) / len(SAMPLES) ))

	# again remove rows with only one genotype because of integrating information from the two sources
	validQTLsnps = np.where([len(set(x))>2 for x in np.array(GT_MERGED_DAT[SAMPLES])])[0]
	GT_MERGED_DAT = GT_MERGED_DAT.iloc[validQTLsnps].reset_index(drop=True)

	#maf = GT_MERGED_DAT[SAMPLES].apply(lambda x: np.sum(x[x!=-1]) / np.sum(x!=-1) / 2, axis = 1)
	#maf = maf.apply(lambda x: 1 - np.max([x, 1-x]))
	#GT_MERGED_DAT = GT_MERGED_DAT.iloc[np.where(maf > 0.05)[0]]

        # not use weights
        WEIGHT_DAT = GT_MERGED_DAT.copy()
        WEIGHT_DAT[SAMPLES] = 1
        compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, GT_MERGED_DAT, WEIGHT_DAT, QTL_dir, saveSuffix = '_withImputation_noWeight')

    else:
        # use weights     
        WEIGHT_DAT = readin_genotype_info(gt_dat=GT_DAT, VCF_dir = VCF_dir, chromosome=CHROMOSOME, samples=SAMPLES)
        compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, GT_DAT, WEIGHT_DAT, QTL_dir)

    ### Compute QTLs using genotype from WGS
    [WGS_DAT, SAMPLES] = read_in_WGS_GT('1k_genome_chr%d' % CHROMOSOME, SAMPLES, WGS_dir, snps = np.array(GT_DAT['CHR_POS']))
    WEIGHT_DAT = WGS_DAT.copy()
    WEIGHT_DAT[SAMPLES] = 1

    compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, WGS_DAT, WEIGHT_DAT, QTL_dir, saveSuffix='_withImputation_realGT')


