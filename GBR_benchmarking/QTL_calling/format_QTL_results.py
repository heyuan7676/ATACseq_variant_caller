import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import sys

def readin_chromosomes(QTL_dir, minDP, cisDist, method):
    try:
        df = pd.read_csv('%s/matrixeQTL_minDP%d_cisDist_%dkb_%s.txt' % (QTL_dir, minDP, (cisDist)/1000, method), sep='\t')
        df_peakLevel = pd.read_csv('%s/matrixeQTL_minDP%d_cisDist_%dkb_%s_peakLevel.txt' % (QTL_dir, minDP, (cisDist)/1000, method), sep='\t')
        return [df, df_peakLevel]
    
    except:
        print('Compute peal-level statistics')
        
    df = pd.DataFrame()
    df_peakLevel = pd.DataFrame()
    for chromosome in range(1,23):
        print('chromsome%d' % chromosome)
        fn1 = '%s/matrixeQTL_chromosome%d_minDP%d_cisDist_%dkb_%s.txt' % (QTL_dir, chromosome, minDP, (cisDist)/1000, method)
        try:
            dfi = pd.read_csv(fn1, sep='\t')
            dfi_peakLevel = dfi.loc[dfi.groupby('gene')["p-value"].idxmin()]
            peak_test_number = pd.DataFrame(dfi.groupby('gene').size())
            peak_test_number.columns = ['Test']
            dfi_peakLevel = peak_test_number.merge(dfi_peakLevel, on = ['gene'])
            dfi_peakLevel['Peak_level_Pvalue'] = [np.min([1,x]) for x in np.array(dfi_peakLevel['p-value']) * np.array(dfi_peakLevel['Test'])]
        except:
            print(fn1)
            print('Result for chr%d not exist' % chromosome)
            continue
            
        df = df.append(dfi)
        df_peakLevel = df_peakLevel.append(dfi_peakLevel)
        
    df['BH_pvalue'] = multipletests(df['p-value'], method='fdr_bh')[1]
    df_peakLevel['BH_pvalue'] = multipletests(df_peakLevel['p-value'], method='fdr_bh')[1]
    
    df = df[['gene', 'SNP', 'p-value', 'BH_pvalue']]
    df.columns = ['gene', 'SNP', 'QTL_level_Pvalue_%s' % method, 'BH_pvalue_%s' % method]
    
    df_peakLevel = df_peakLevel[['gene', 'SNP', 'p-value','Peak_level_Pvalue', 'BH_pvalue']]
    df_peakLevel.columns = ['gene', 'SNP', 'QTL_level_Pvalue_%s' % method,'Peak_level_Pvalue_%s' % method, 'BH_pvalue_%s' % method]
    
    df.to_csv('%s/matrixeQTL_minDP%d_cisDist_%dkb_%s.txt' % (QTL_dir, minDP, (cisDist)/1000, method), sep='\t', index=False)
    df_peakLevel.to_csv('%s/matrixeQTL_minDP%d_cisDist_%dkb_%s_peakLevel.txt' % (QTL_dir, minDP, (cisDist)/1000, method), sep='\t', index=False)
    

    return [df, df_peakLevel]



minDP = 3
QTL_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/QTLs_MACS2'


for cisDist in [0, 1000, 10000, 100000, 1000000]:
    for method in ['WGS', 'GC', 'Imputation', 'Integration']:
        print(method)
        readin_chromosomes(QTL_dir, minDP, cisDist, method)






