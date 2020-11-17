import sys
import time
import os
import pdb
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from statsmodels.regression import linear_model as sm
from prepare_matrixQTL import *
from QTL_calling import *


def read_in_peaks_Zenodo(PEAK_dir, chromosome):
    peak_dat = pd.read_csv('%s/chromosome%d_corrected_fpkm.txt' % (PEAK_dir, chromosome), sep='\t')
    peak_dat.columns = ['PEAK'] + list(peak_dat.columns[1:])

    peak_loc = pd.read_csv('%s/chromosome%d_loc.bed' % (PEAK_dir, chromosome), sep='\t')
    peak_loc.columns = ['PEAK', 'CHR', 'START', 'END']

    peak_dat = peak_dat.merge(peak_loc, on = 'PEAK')

    samples = [x for x in peak_dat.columns if x.startswith('HG')]
    return [peak_dat, samples]



if __name__ == "__main__":
    CHROMOSOME = int(sys.argv[1])
    WINDOW = int(sys.argv[2])
    SUFFIX = '_minDP2'

    PEAK_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks/Zenodo'

    # root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_subsample_0.5'
    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie'
    Genotype_dir = '%s/Called_GT/' % root_dir
    VCF_dir = '%s/VCF_files' % root_dir
    QTL_dir = '%s/QTLs/Zenodo' % root_dir
    #os.makedirs(QTL_dir, exist_ok = True)

    print('Call ac-QTLs for chromosome %d with window = %dkb' % (CHROMOSOME, WINDOW/1000))

    ## read in data
    [PEAK_DAT, SAMPLES] = read_in_peaks_Zenodo(PEAK_dir, CHROMOSOME)

    ## read in data
    GT_DAT = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix=SUFFIX)
    #WEIGHT_DAT = readin_genotype_info(gt_dat=GT_DAT, VCF_dir = VCF_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix=SUFFIX)

    ## compute QTLs
    #compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, GT_DAT, WEIGHT_DAT, QTL_dir, suffix=SUFFIX)

    WEIGHT_DAT = GT_DAT.copy()
    WEIGHT_DAT[SAMPLES] = 1

    ## compute QTLs
    compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, GT_DAT, WEIGHT_DAT, QTL_dir, suffix=SUFFIX, saveSuffix = '_noWeight')

