import numpy as np
import pandas as pd
import time

import sys
import os
import pdb

from prepare_matrixQTL import *
from QTL_calling import *

def read_in_peaks(PEAK_dir, chromosome):
    peak_dat = pd.read_csv('%s/peak_by_sample_matrix_RPKM_corrected_chromosome%d.txt' % (PEAK_dir, chromosome), sep='\t')
    samples = [x for x in peak_dat.columns if x.startswith('HG')]

    return [peak_dat, samples]


def read_in_peaks_Zenodo(PEAK_dir, chromosome):
    peak_dat = pd.read_csv('%s/chromosome%d_corrected_fpkm.txt' % (PEAK_dir, chromosome), sep='\t')
    peak_dat.columns = ['PEAK'] + list(peak_dat.columns[1:])

    peak_loc = pd.read_csv('%s/chromosome%d_loc.bed' % (PEAK_dir, chromosome), sep='\t')
    peak_loc.columns = ['PEAK', 'CHR', 'START', 'END']

    peak_dat = peak_dat.merge(peak_loc, on = 'PEAK')

    samples = [x for x in peak_dat.columns if x.startswith('HG')]
    return [peak_dat, samples]



def read_in_WGS_GT(prefix, samples_peaks, use_only_called_genotypes = False):
    WGS_fn = '%s/%s.genotypes.tsv' % (WGS_dir, prefix)
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t')
    WGS_result = WGS_result.drop_duplicates()
    
    samples = [x for x in WGS_result.columns if x.startswith('HG')]
    WGS_result.index = WGS_result[['#CHROM', 'POS']].apply(lambda x: '_'.join((str(x[0]), str(x[1]))), axis=1)

    WGS_result = WGS_result.replace('./.', '0')
    WGS_result = WGS_result.replace(0, '0')

    # remove rows with less than 3 samples
    valid_snps = np.where(np.sum(np.array(WGS_result[samples]) != '0', axis=1) >= 3)[0]
    WGS_result = WGS_result.iloc[valid_snps].reset_index(drop=True)

    # remove rows with only one genotype
    validQTLsnps = np.where([len(set(x))>2 for x in np.array(WGS_result[samples])])[0]
    WGS_result = WGS_result.iloc[validQTLsnps].reset_index(drop=True)

    if use_only_called_genotypes: 
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
    SUFFIX = '_minDP2'
    #SUFFIX = '_minDP2_noWeight'

    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/WithInPeaks_Bowtie'
    PEAK_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks'

    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie'
    QTL_dir = '%s/QTLs' % root_dir
    Genotype_dir = '%s/Called_GT/' % root_dir

    [PEAK_DAT, samples_peaks] = read_in_peaks(PEAK_dir, CHROMOSOME)
   
    #[WGS_dat, SAMPLES] = read_in_WGS_GT('1k_genome_chr%d' % CHROMOSOME, samples_peaks, use_only_called_genotypes = True)
    [WGS_dat, SAMPLES] = read_in_WGS_GT('chr%d.WithInPeaks' % CHROMOSOME, samples_peaks)
    WGS_DAT = obtain_numerical_gt(WGS_dat, SAMPLES)
    # again remove rows with only one genotype (ie. A/T and T/A)
    validQTLsnps = np.where([len(set(x))>2 for x in np.array(WGS_DAT[SAMPLES])])[0]
    WGS_DAT = WGS_DAT.iloc[validQTLsnps].reset_index(drop=True)

    WEIGHT_DAT = WGS_DAT.copy()
    WEIGHT_DAT[SAMPLES] = 1
  
    save_dir = os.path.join(QTL_dir, 'WithInPeaks_Bowtie')
    try:
        os.makedirs(save_dir) 
    except:
        pass
    compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, WGS_DAT, WEIGHT_DAT, save_dir=save_dir, suffix=SUFFIX, saveSuffix='_realGT')
    

