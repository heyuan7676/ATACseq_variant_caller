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
    
    SNPs_close = np.where((np.array(gt_numerical_dat['POS']) < int(end + WINDOW)) * (np.array(gt_numerical_dat['POS']) > int(start - WINDOW)))[0]
   
    samples = [x for x in gt_numerical_dat.columns if x.startswith('HG')] 
    bb = np.array(gt_numerical_dat.iloc[SNPs_close][samples])
    cc = np.array(post_pp_dat.iloc[SNPs_close][samples])
    aa = np.repeat(np.array(peaki[samples])[np.newaxis], len(bb), axis=0)

    datapoints = zip(aa,bb,cc)
    try:
    	pvalues = list(map(compute_QTL_gti_peaki, datapoints)) 
    except:
	pdb.set_trace()
    QTL_result_peaki = pd.DataFrame({"PeakID": peaki['PEAK'], "Start": peaki['START'], "End": peaki['END'], "CHR_POS": gt_numerical_dat.iloc[SNPs_close]['CHR_POS'], "POS": gt_numerical_dat.iloc[SNPs_close]['POS'], "P-value": pvalues})

    return QTL_result_peaki



def compute_QTLs(chromosome, WINDOW, peak_df, genotype_df, weight_df, save_dir, saveSuffix = ''):

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
    QTL_results.to_csv('%s/CHR%d_caQTLs_WINDOW_%skb%s.txt' % (save_dir, chromosome, str(WINDOW/1000.0), saveSuffix), sep='\t', index = False)


