import numpy as np
import pandas as pd

import sys
import os

from prepare_data_matrix import *
from QTL_calling import *
from statsmodels.stats.multitest import multipletests
from Derive_QTLs import *


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
        captured_QTLs = readin_QTL_results(QTL_dir, VCF_dir, SAMPLEs, WINDOW, suffix = suffix)
        obtain_performance(captured_QTLs, alpha, '%s/%s_%s_%fkb%s.txt' % (save_dir, peak_calling, GT_subDir, str(WINDOW/1000.0), suffix))




