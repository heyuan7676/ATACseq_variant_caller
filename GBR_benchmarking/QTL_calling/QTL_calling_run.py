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
    peak_calling = sys.argv[3]
    WINDOW = int(sys.argv[4])
    minDP = int(sys.argv[5])
    method = sys.argv[6]

    alignment_dir = 'alignment_bowtie' # alignment_subsample_0.5
    GT_subDir = 'minDP%d' % minDP

    root_dir = '%s/%s' % (ROOT_DIR, alignment_dir)
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'

    if peak_calling == 'macs2':
    	PEAK_dir = '%s/Peaks_MACS2' % root_dir
        QTL_dir = '%s/QTLs_MACS2/%s' % (root_dir, GT_subDir)

    if peak_calling == 'macs2_combined':
        PEAK_dir = '%s/Peak_version1/Peaks/combined' % root_dir
        QTL_dir = '%s/QTLs_MACS2_combined/%s' % (root_dir, GT_subDir)

    elif peak_calling == 'Genrich':
        PEAK_dir = '%s/Peaks_Genrich' % root_dir
        QTL_dir = '%s/QTLs_Genrich/%s' % (root_dir, GT_subDir)

    elif peak_calling == 'Genrich_combined':
        PEAK_dir = '%s/Peaks_Genrich/combined' % root_dir
        QTL_dir = '%s/QTLs_Genrich_combined/%s' % (root_dir, GT_subDir)

    if not os.path.exists(QTL_dir):
        os.makedirs(QTL_dir)

    QTL_dir_imputed = '%s/Imputation' % QTL_dir
    if not os.path.exists(QTL_dir_imputed):
	os.makedirs(QTL_dir_imputed)

    QTL_dir_integration = '%s/Integration' % QTL_dir
    if not os.path.exists(QTL_dir_integration):
        os.makedirs(QTL_dir_integration)

    ## read in data
    [PEAK_DAT, SAMPLES_Peaks] = read_in_peaks(PEAK_dir, CHROMOSOME)

    ## align the samples with samples from WGS
    samples = pd.read_csv('../Test_samples.txt', sep='\t', header = None)
    samples = np.array(samples[0])

    WGS_fn = '%s/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.GT.FORMAT' % (WGS_dir, 22)
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t', low_memory=False, nrows = 10)
    SAMPLES_WGS  = [x for x in WGS_result.columns if x.startswith('HG')]
    SAMPLES = list(np.intersect1d(SAMPLES_Peaks, SAMPLES_WGS))
    SAMPLES = list(np.intersect1d(SAMPLES, samples))
    print('%d samples are used for QTL analysis' % len(SAMPLES))
    print('Call ca-QTLs for chromosome %d with window = %skb; using peaks from %s, variants from %s\n' % (CHROMOSOME, str(WINDOW/1000.0), peak_calling, method))

    flag = 0

    if method == 'Integration':
        fn = os.path.join(QTL_dir_integration, 'CHR%d_caQTLs_WINDOW_%skb%s.txt' % (CHROMOSOME, str(WINDOW/1000.0), '_Integration'))
        try:
            dat = pd.read_csv(fn, nrows=10, sep='\t')
            flag = 1
        except:
            print('Call ca-QTLs for chromosome %d with window = %skb' % (CHROMOSOME, str(WINDOW/1000.0)))

        if flag:
            sys.exit()
	else:
	    ####### use both
    	    print('Read in data from both sources...')
            Genotype_dir = '%s/Integration/%s' % (root_dir, GT_subDir)
            [GT_DAT_Integrated, numbers_integrated] = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix = '.with_Inconsistent')

            compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, GT_DAT_Integrated, QTL_dir_integration, saveSuffix = '_Integration')



    if method == 'Imputation':
	fn = os.path.join(QTL_dir_imputed, 'CHR%d_caQTLs_WINDOW_%skb%s.txt' % (CHROMOSOME, str(WINDOW/1000.0), '_Imputation'))
        try:
            dat = pd.read_csv(fn, nrows=10, sep='\t')
            flag = 1
        except:
            print('Call ca-QTLs for chromosome %d with window = %skb' % (CHROMOSOME, str(WINDOW/1000.0)))

        if flag:
            sys.exit()
        else:
	    ####### use imputation
    	    print('Read in genotype data from imputation...')
    	    Genotype_dir = '%s/Imputation/%s' % (root_dir, GT_subDir)
    	    [GT_DAT_Imputed, numbers_imputation] = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix = '')

	    compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, GT_DAT_Imputed, QTL_dir_imputed, saveSuffix = '_Imputation')



    if method == 'Called':
	fn = os.path.join(QTL_dir, 'CHR%d_caQTLs_WINDOW_%skb%s.txt' % (CHROMOSOME, str(WINDOW/1000.0), '_Called'))
        try:
            dat = pd.read_csv(fn, nrows=10, sep='\t')
            flag = 1
        except:
            print('Call ca-QTLs for chromosome %d with window = %skb' % (CHROMOSOME, str(WINDOW/1000.0)))

        if flag:
            sys.exit()
        else:
	    ## read in data
	    print('Read in genotype data from ATAC-seq reads...')
	    Genotype_dir = '%s/VCF_files/%s' % (root_dir, GT_subDir)
	    [GT_DAT, numbers_atac] = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix = '.recode')
    	    compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, GT_DAT, QTL_dir, saveSuffix = '_Called')


    if 0:
        numbers = [numbers_atac, numbers_imputation]
        numbers = pd.DataFrame(numbers)
        numbers.columns = ["All_variants_called", "Variants_with_more_than_one_genotype", "Variants_with_MinorAC_>=_3", "Variants_with_MajoyAC_>=_3"]
        numbers.index = ['ATAC_reads', 'imputation']
        numbers.to_csv('variants_numbers/CHR%d_WINDOW_%skb_%s_variants_Number.txt' % (CHROMOSOME, str(WINDOW/1000.0), peak_calling), sep='\t')


    save_matrix = False
    if save_matrix:
        save_dir = os.path.join(QTL_dir, 'Save_Matrix')
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        GT_DAT[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d.txt' % (save_dir, CHROMOSOME), sep='\t', index = False)
        GT_DAT[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d_loc.bed' % (save_dir, CHROMOSOME), sep='\t', index=False)


        save_dir = os.path.join(QTL_dir_imputed, 'Save_Matrix')
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        GT_DAT_Imputed[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d.txt' % (save_dir, CHROMOSOME), sep='\t', index = False)
        GT_DAT_Imputed[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d_loc.bed' % (save_dir, CHROMOSOME), sep='\t', index=False)


        save_dir = os.path.join(QTL_dir_integration, 'Save_Matrix')
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        GT_DAT_Integrated[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d.txt' % (save_dir, CHROMOSOME), sep='\t', index = False)
        GT_DAT_Integrated[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d_loc.bed' % (save_dir, CHROMOSOME), sep='\t', index=False)
#
