import sys
import time
import os
import pdb
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from statsmodels.regression import linear_model as sm
from scipy import stats

#### Compute QTL

# Genotype data: gt_numerical_dat
# Peak data: RPKM_dat
# Genotype data weights: post_pp_dat


def compute_QTL_gti_peaki(datapoint):
    #gt_sample = np.array(gt_numerical_dat[samples].iloc[0])
    #peak_sample = np.array(RPKM_dat[samples].iloc[0])
    #weight_sample = np.array(post_pp_dat[samples].iloc[0])

    [peak_sample, gt_sample, weight_sample] = datapoint
    valid_samples = np.where(gt_sample!= -1)[0]
    
    y = np.array(peak_sample[valid_samples])
    y = y.astype(float)
    x = np.array(gt_sample[valid_samples])
    x_weights = np.array(weight_sample[valid_samples])
   
    x = sm.add_constant(x) 
    wls_model = sm.WLS(y, x, weights = x_weights)
    results = wls_model.fit()

    return results.pvalues[1]



def compute_QTL_peaki(WINDOW, peaki, gt_numerical_dat, post_pp_dat):
    QTL_result_peaki = []
    [start, end] = [peaki['START'], peaki['END']]
    
    SNPs_close = (np.array(gt_numerical_dat['POS']) < end + WINDOW) & (np.array(gt_numerical_dat['POS']) > start - WINDOW)
    SNPs_close = np.where(SNPs_close)[0]
   
    samples = [x for x in gt_numerical_dat.columns if x.startswith('HG')] 
    bb = np.array(gt_numerical_dat.iloc[SNPs_close][samples])
    cc = np.array(post_pp_dat.iloc[SNPs_close][samples])
    aa = np.repeat(np.array(peaki[samples])[np.newaxis], len(bb), axis=0)

    datapoints = zip(aa,bb,cc)
    pvalues = list(map(compute_QTL_gti_peaki, datapoints)) 
    QTL_result_peaki = pd.DataFrame({"PeakID": peaki['PEAK'], "CHR_POS": gt_numerical_dat.iloc[SNPs_close]['CHR_POS'], "P-value": pvalues})

    return QTL_result_peaki



def compute_QTLs(chromosome, WINDOW, peak_df, genotype_df, weight_df, save_dir, suffix = '', saveSuffix = ''):

    QTL_results = pd.DataFrame()
    print('Compute QTLs for %d peaks ...' % len(peak_df))
    start = time.time()
    for p in range(len(peak_df)):
        token = compute_QTL_peaki(WINDOW, peak_df.iloc[p], genotype_df, weight_df)
        QTL_results = QTL_results.append(token)
        if (p+1) % 500 == 0:
            print('    %d peaks finished' % (p+1))
    end = time.time()
    print('Tested %d associations for chromosome%d' % (len(QTL_results), chromosome))
    print('    Used %f miniutes' % ((end - start)/60))
    print("") 
    QTL_results = QTL_results.sort_values('P-value')
    QTL_results.to_csv('%s/CHR%d_acQTLs_WINDOW_%dkb%s%s.txt' % (save_dir, chromosome, WINDOW/1000, suffix, saveSuffix), sep='\t', index = False)


