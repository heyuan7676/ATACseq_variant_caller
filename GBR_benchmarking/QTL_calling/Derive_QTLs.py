import numpy as np
import pandas as pd

import sys
import os

from prepare_data_matrix import *
from QTL_calling import *
from statsmodels.stats.multitest import multipletests



def readin_both_results(QTL_dir, chromsome, window, suffix):
    fn1 = '%s/CHR%d_caQTLs_WINDOW_%skb%s.txt' % (QTL_dir, chromsome, str(window/1000.0), '_realGT')
    real_QTL = pd.read_csv(fn1, sep='\t')
    
    fn2 = '%s/CHR%d_caQTLs_WINDOW_%skb%s.txt' % (QTL_dir, chromsome, str(window/1000.0), suffix)
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





def readin_QTL_results(QTL_dir, VCF_dir, SAMPLEs, window, suffix):

    QTLs_dat = pd.DataFrame()
    QTLs_dat_PP = pd.DataFrame()

    for CHROMOSOME in range(1,22):
        try:
            token = readin_both_results(QTL_dir, CHROMOSOME, window, suffix)
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



def Derive_QTLs(df):

    group = df.copy()
    
    QTLs_real = group.loc[group.groupby('PeakID').pvalue_real_Peak.idxmin()]
    QTLs_real['BH_pvalue_real'] = multipletests(QTLs_real['pvalue_real_Peak'], method='fdr_bh')[1]

    QTLs_called = group.loc[group.groupby('PeakID').pvalue_called_Peak.idxmin()]
    QTLs_called['BH_pvalue_called'] = multipletests(QTLs_called['pvalue_called_Peak'], method='fdr_bh')[1]

    return [QTLs_real, QTLs_called]


