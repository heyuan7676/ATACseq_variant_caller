import sys
import time
import os
import pdb
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from statsmodels.regression import linear_model as sm
from prepare_data_matrix import *
from QTL_calling import *

def read_in_peaks(PEAK_dir, chromosome):
    peak_dat = pd.read_csv('%s/peak_by_sample_matrix_RPKM_corrected_chromosome%d.txt' % (PEAK_dir, chromosome), sep='\t')
    samples = [x for x in peak_dat.columns if x.startswith('HG')]
    return [peak_dat, samples]


if __name__ == "__main__":
    CHROMOSOME = int(sys.argv[1])
    WINDOW = int(sys.argv[2])
    useWeight = sys.argv[3]

    alignment_dir = 'alignment_bowtie' # alignment_subsample_0.5
    peak_calling = 'macs2'
    GT_subDir = 'minDP2'
    SUFFIX = '_%s' % GT_subDir

    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/%s' % alignment_dir
    VCF_dir = '%s/VCF_files' % root_dir
    Genotype_dir = '%s/Called_GT/%s' % (root_dir, GT_subDir)

    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'

    if peak_calling == 'macs2':
    	PEAK_dir = '%s/Peaks' % root_dir
        QTL_dir = '%s/QTLs/%s' % (root_dir, GT_subDir)
        os.makedirs(QTL_dir)

    elif peak_calling == 'Genrich':
        PEAK_dir = '%s/Peaks_Genrich' % root_dir
        QTL_dir = '%s/QTLs_Genrich/%s' % (root_dir, GT_subDir)
        os.makedirs(QTL_dir)


    print('Call ca-QTLs for chromosome %d with window = %dkb' % (CHROMOSOME, WINDOW/1000))

    ## read in data
    [PEAK_DAT, SAMPLES_Peaks] = read_in_peaks(PEAK_dir, CHROMOSOME)

    ## align the samples with samples from WGS
    WGS_fn = '%s/%s.genotypes.tsv' % (WGS_dir, '1k_genome_chr22')
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t', nrows = 10)
    SAMPLES_WGS = [x for x in WGS_result.columns if x.startswith('HG')]
    
    SAMPLES = list(np.intersect1d(SAMPLES_Peaks, SAMPLES_WGS))

    print('%d samples are used for QTL analysis' % len(SAMPLES))

    save_dir = os.path.join(QTL_dir)
    try:
        os.makedirs(save_dir)
    except:
        pass

    ## read in data
    GT_DAT = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix=SUFFIX)
    if useWeight:
        WEIGHT_DAT = readin_genotype_info(gt_dat=GT_DAT, VCF_dir = VCF_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix=SUFFIX)
        ## compute QTLs
        compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, GT_DAT, WEIGHT_DAT, QTL_dir, suffix=SUFFIX)

    else:
        WEIGHT_DAT = GT_DAT.copy()
        WEIGHT_DAT[SAMPLES] = 1
        ## compute QTLs
        compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, GT_DAT, WEIGHT_DAT, QTL_dir, suffix=SUFFIX, saveSuffix = '_noWeight')

    ### Compute QTLs using genotype from WGS
    [WGS_DAT, SAMPLES] = read_in_WGS_GT('1k_genome_chr%d' % CHROMOSOME, samples_peaks, WGS_dir, snps = np.array(GT_DAT['CHR_POS']))
    WEIGHT_DAT = WGS_DAT.copy()
    WEIGHT_DAT[SAMPLES] = 1

    compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, WGS_DAT, WEIGHT_DAT, save_dir=save_dir, suffix=SUFFIX, saveSuffix='_realGT')


