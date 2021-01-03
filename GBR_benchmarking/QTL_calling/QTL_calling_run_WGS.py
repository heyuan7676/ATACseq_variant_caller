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

    alignment_dir = 'alignment_bowtie' # alignment_subsample_0.5
    GT_subDir = 'minDP2'

    root_dir = '%s/%s' % (ROOT_DIR, alignment_dir)
    VCF_dir = '%s/VCF_files' % root_dir
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

    fn = os.path.join(QTL_dir, 'CHR%d_caQTLs_WINDOW_%skb%s.txt' % (CHROMOSOME, str(WINDOW/1000.0), '_realGT_all'))
    flag = 0
    try:
	dat = pd.read_csv(fn, nrows=10, sep='\t')
	flag = 1
    except:
	print(fn)
        print('Call ca-QTLs for chromosome %d with window = %skb' % (CHROMOSOME, str(WINDOW/1000.0)))

    if flag:
	sys.exit()

    ## use only the Test samples
    samples = pd.read_csv('../Test_samples.txt', sep='\t', header = None)
    samples = np.array(samples[0])

    ## read in data
    [PEAK_DAT, SAMPLES_Peaks] = read_in_peaks(PEAK_dir, CHROMOSOME)
    samples = np.intersect1d(samples, SAMPLES_Peaks)

    ### Compute QTLs using genotype from WGS
    [WGS_DAT, SAMPLES, numbers_WGS_called] = read_in_WGS_GT(WGS_dir, CHROMOSOME, samples)
    compute_QTLs(CHROMOSOME, WINDOW, PEAK_DAT, WGS_DAT, QTL_dir, saveSuffix='_realGT_all')

    numbers = [numbers_WGS_called]
    numbers = pd.DataFrame(numbers)
    numbers.columns = ["All_variants_called", "Variants_with_more_than_one_genotype", "Variants_with_MAC_>=_3", "Variants_with_MajoyAC_>=_3"]
    numbers.index = ['genotype']
    numbers.to_csv('variants_numbers/CHR%d_WINDOW_%skb_%s_variants_Number_WGS.txt' % (CHROMOSOME, str(WINDOW/1000.0), peak_calling), sep='\t')



