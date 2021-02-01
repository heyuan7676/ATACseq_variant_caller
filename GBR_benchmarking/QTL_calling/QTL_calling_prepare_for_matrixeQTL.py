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


if __name__ == "__main__":
    minDP = int(sys.argv[1])
    CHROMOSOME = int(sys.argv[2])
    GT_subDir = 'minDP%d' % minDP
    
    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie'
    save_GT_dir = '%s/Dosage_for_QTL/%s' % (root_dir, GT_subDir)

    GT_dir = os.path.join(save_GT_dir, 'GC')
    GT_dir_imputed = os.path.join(save_GT_dir, 'Imputation')
    GT_dir_integration = os.path.join(save_GT_dir, 'Integration')
    os.makedirs(GT_dir, exist_ok = True)
    os.makedirs(GT_dir_imputed, exist_ok = True)
    os.makedirs(GT_dir_integration, exist_ok = True)

    os.makedirs(os.path.join(GT_dir, 'Save_Matrix'), exist_ok = True)
    os.makedirs(os.path.join(GT_dir_imputed, 'Save_Matrix'), exist_ok = True)
    os.makedirs(os.path.join(GT_dir_integration, 'Save_Matrix'), exist_ok = True)

    ## align the samples with samples from WGS
    samples = pd.read_csv('../samples.txt', sep='\t', header = None)
    SAMPLES = list(np.sort(np.array(samples[0])))
    print('%d samples are used for QTL analysis' % len(SAMPLES))

    ####### use both
    print('Read in data from both sources...')
    Genotype_dir = '%s/Integration/%s' % (root_dir, GT_subDir)
    [GT_DAT_Integrated, numbers_integrated] = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix = '.THR0.0_Y_random_forest')
    save_dir = os.path.join(GT_dir_integration, 'Save_Matrix')
    GT_DAT_Integrated[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d.txt' % (save_dir, CHROMOSOME), sep='\t', index = False)
    GT_DAT_Integrated[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d_loc.bed' % (save_dir, CHROMOSOME), sep='\t', index=False)


    ####### use imputation
    print('Read in genotype data from imputation...')
    Genotype_dir = '%s/Imputation/%s' % (root_dir, GT_subDir)
    [GT_DAT_Imputed, numbers_imputation] = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix = '')
    save_dir = os.path.join(GT_dir_imputed, 'Save_Matrix')
    GT_DAT_Imputed[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d.txt' % (save_dir, CHROMOSOME), sep='\t', index = False)
    GT_DAT_Imputed[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d_loc.bed' % (save_dir, CHROMOSOME), sep='\t', index=False)


    ## read in data
    print('Read in genotype data from ATAC-seq reads...')
    Genotype_dir = '%s/VCF_files/%s' % (root_dir, GT_subDir)
    [GT_DAT, numbers_atac] = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix = '.recode')

    save_dir = os.path.join(GT_dir, 'Save_Matrix')
    GT_DAT[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d.txt' % (save_dir, CHROMOSOME), sep='\t', index = False)
    GT_DAT[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d_loc.bed' % (save_dir, CHROMOSOME), sep='\t', index=False)

