import numpy as np
import pandas as pd

import sys
import os

from prepare_data_matrix import *
from QTL_calling import *
from statsmodels.stats.multitest import multipletests



def readin_both_results(chromsome, window, suffix):
    fn1 = '%s/CHR%d_caQTLs_WINDOW_%dkb%s.txt' % (QTL_dir, chromsome, window/1000, '_realGT')
    real_QTL = pd.read_csv(fn1, sep='\t')
    
    fn2 = '%s/CHR%d_caQTLs_WINDOW_%dkb%s.txt' % (QTL_dir, chromsome, window/1000, suffix)
    call_QTL = pd.read_csv(fn2, sep='\t')

    real_QTL.columns = ['CHR_POS','pvalue_real', 'PeakID']
    call_QTL.columns = ['CHR_POS','pvalue_called', 'PeakID']

    QTLs_in_both = real_QTL.merge(call_QTL, on = ['CHR_POS', 'PeakID'])
    
    return QTLs_in_both




def readin_both_results_matrixeQTL(chromsome, window, suffix):
    fn1 = '%s/matrixeQTL_chromosome%d%s_cisDist%dkb%s.txt' % (QTL_dir, chromsome, suffix, window/1000, '_realGT')
    real_QTL = pd.read_csv(fn1, sep='\t', usecols=[0, 1, 4])
    
    fn2 = '%s/matrixeQTL_chromosome%d%s_cisDist%dkb%s.txt' % (QTL_dir, chromsome, suffix, window/1000, '')
    call_QTL = pd.read_csv(fn2, sep='\t', usecols=[0, 1, 4])

    real_QTL.columns = ['CHR_POS', 'PeakID','pvalue_real']
    call_QTL.columns = ['CHR_POS', 'PeakID','pvalue_called']

    QTLs_in_both = real_QTL.merge(call_QTL, on = ['CHR_POS', 'PeakID'])
    
    return QTLs_in_both





def readin_QTL_results(SAMPLEs, suffix):

    QTLs_dat = pd.DataFrame()
    QTLs_dat_PP = pd.DataFrame()

    for CHROMOSOME in range(1,22):
        try:
            token = readin_both_results(CHROMOSOME, WINDOW, suffix)
            QTLs_dat = QTLs_dat.append(token)
        except:
            print('chromosome%d not exist' % CHROMOSOME)
            continue
   
        token_pp = readin_genotype_info(token, VCF_dir, CHROMOSOME, SAMPLEs)
        QTLs_dat_PP = QTLs_dat_PP.append(token_pp[SAMPLEs])

    QTLs_dat['min_nonzero_PP'] = [np.min(x[x>0]) for x in np.array(QTLs_dat_PP)]

    QTLs_dat['PP_category'] = '0'
    QTLs_dat.loc[QTLs_dat['min_nonzero_PP'] < 0.4, 'PP_category'] = '1'
    QTLs_dat.loc[(QTLs_dat['min_nonzero_PP'] < 0.5) & (QTLs_dat['min_nonzero_PP'] >= 0.4), 'PP_category'] = '2'
    QTLs_dat.loc[(QTLs_dat['min_nonzero_PP'] < 1) & (QTLs_dat['min_nonzero_PP'] >= 0.5), 'PP_category'] = '3'
    QTLs_dat.loc[QTLs_dat['min_nonzero_PP'] == 1, 'PP_category'] = '4'

    number_test_per_peak = pd.DataFrame(QTLs_dat.groupby('PeakID').size())
    number_test_per_peak.columns = ['Number_tests']
    QTLs_dat = QTLs_dat.merge(number_test_per_peak, left_on='PeakID', right_index=True)    
    QTLs_dat['pvalue_real_Peak'] = list(map(lambda x: min(1, x), list(QTLs_dat['pvalue_real'] * QTLs_dat['Number_tests'])))
    QTLs_dat['pvalue_called_Peak'] = list(map(lambda x: min(1, x), list(QTLs_dat['pvalue_called'] * QTLs_dat['Number_tests'])))

    return QTLs_dat



def obtain_performance(QTLs_dat, alpha, save_file):
    collect_metric = []

    QTLs_dat['PP_category'] = [int(x) for x in QTLs_dat['PP_category']]
    color_group = {"1": ["0<minPP<0.4", 'red'],
                   "2": ["0.4=<minPP<0.5", 'orange'],
                   "3": ["0.5<=minPP<1", 'blue'],
                   "4": ["minPP=1", 'green']}

    for thr in [1,2,3,4]:
        group = QTLs_dat[QTLs_dat['PP_category'] >= thr]
    
        QTLs_real = group.loc[group.groupby('PeakID').pvalue_real_Peak.idxmin()]
        A = QTLs_real.iloc[np.where(multipletests(QTLs_real['pvalue_real_Peak'], method='fdr_bh', alpha=alpha)[0])[0]]
    
        QTLs_called = group.loc[group.groupby('PeakID').pvalue_called_Peak.idxmin()]
        B = QTLs_called.iloc[np.where(multipletests(QTLs_called['pvalue_called_Peak'], method='fdr_bh', alpha=alpha)[0])[0]]
    
        A = np.array(A.apply(lambda x: '%s:%s' % ((x['CHR_POS'], x['PeakID'])), axis=1))
        B = np.array(B.apply(lambda x: '%s:%s' % ((x['CHR_POS'], x['PeakID'])), axis=1))
   
        collect_metric.append([color_group[str(thr)][0],
                               len(set(A)), 
                               len(set(B)), 
                               len(set(np.intersect1d(A, B)))])
    
    collect_metric_df = pd.DataFrame(collect_metric)
    collect_metric_df.columns = ["GQ_threshold","True_QTL", "Called_QTL", "Overlap"]
    collect_metric_df['Recall'] = collect_metric_df['Overlap'] / collect_metric_df['True_QTL']
    collect_metric_df['Precision'] = collect_metric_df['Overlap'] / collect_metric_df['Called_QTL']
    collect_metric_df['GQ_threshold'] = ['PP>0', 'PP>0.4', 'PP>0.5', 'PP=1']
    collect_metric_df.to_csv(save_file, sep='\t', index = False)


if __name__ == '__main__':
    peak_calling = sys.argv[1]
    WINDOW = 0

    alignment_dir = 'alignment_bowtie' # alignment_subsample_0.5
    GT_subDir = 'minDP2'

    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/%s' % alignment_dir
    VCF_dir = '%s/VCF_files' % root_dir
    Genotype_dir = '%s/Called_GT/%s' % (root_dir, GT_subDir)

    if peak_calling == 'macs2':
        PEAK_dir = '%s/Peaks' % root_dir
        QTL_dir = '%s/QTLs/%s' % (root_dir, GT_subDir)

    if peak_calling == 'macs2_combined':
        PEAK_dir = '%s/Peaks/combined' % root_dir
        QTL_dir = '%s/QTLs_combined/%s' % (root_dir, GT_subDir)

    elif peak_calling == 'Genrich':
        PEAK_dir = '%s/Peaks_Genrich' % root_dir
        QTL_dir = '%s/QTLs_Genrich/%s' % (root_dir, GT_subDir)

    if not os.path.exists(QTL_dir):
        os.makedirs(QTL_dir)

    save_dir = 'Evaluation_metrics'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    ## read in data
    SAMPLEs_Peaks = readin_peak_samples(PEAK_dir)

    ## align the samples with samples from WGS
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'
    WGS_fn = '%s/%s.genotypes.tsv' % (WGS_dir, '1k_genome_chr22')
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t', nrows = 10)
    SAMPLEs_WGS = [x for x in WGS_result.columns if x.startswith('HG')]

    SAMPLEs = list(np.intersect1d(SAMPLEs_Peaks, SAMPLEs_WGS))

    
    alpha = 0.05
    for suffix in ['', '_noWeight']:
        captured_QTLs = readin_QTL_results(SAMPLEs, suffix = '')
        obtain_performance(captured_QTLs, alpha, '%s/%s_%s%s.txt' % (save_dir, peak_calling, GT_subDir, suffix))




