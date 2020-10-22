import numpy as np
import pandas as pd
import time

import sys
import os

from QTL_calling import *


def read_in_WGS_GT(samples_peaks):

    WGS_fn = '%s/1k_genome_chr%d.genotypes.tsv' % (WGS_dir, CHROMOSOME)
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t')
    WGS_result = WGS_result.drop_duplicates()
    
    samples = [x for x in WGS_result.columns if x.startswith('HG')]
    WGS_result.index = WGS_result[['#CHROM', 'POS']].apply(lambda x: '_'.join((str(x[0]), str(x[1]))), axis=1)
    
    # use only the SNPs from called genotypes
    QTL_by_called_genotypes = pd.read_csv('%s/CHR%d_acQTLs_WINDOW_%dkb%s.txt' % (QTL_dir, CHROMOSOME, WINDOW/1000, SUFFIX), sep='\t')    
    WGS_result = WGS_result.loc[QTL_by_called_genotypes['CHR_POS']]
    
    WGS_result = WGS_result.drop_duplicates()
    WGS_result.columns = ['CHR', 'POS'] + list(WGS_result.columns[2:])
    WGS_result.index.name = 'CHR_POS'
    WGS_result = WGS_result.reset_index()
    
    samples = list(np.intersect1d(samples_peaks, WGS_result.columns))
    WGS_result = WGS_result[list(WGS_result.columns[:3]) + samples]
    
    return [WGS_result, samples]



if __name__ == '__main__':
    CHROMOSOME = int(sys.argv[1])
    WINDOW = int(sys.argv[2])
    SUFFIX = '_GQ1'

    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'
    PEAK_DIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks'
    QTL_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/QTLs'
    Genotype_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Called_GT/'

    [PEAK_DAT, samples_peaks] = read_in_peaks(PEAK_DIR, CHROMOSOME)
    
    [WGS_dat, SAMPLES] = read_in_WGS_GT(samples_peaks)
    WGS_DAT = obtain_numerical_gt(WGS_dat, SAMPLES)
    WEIGHT_DAT = WGS_DAT.copy()
    WEIGHT_DAT[SAMPLES] = 1
    
    compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, WGS_DAT, WEIGHT_DAT, QTL_dir=QTL_dir, suffix=SUFFIX, saveSuffix='_realGT')
    

