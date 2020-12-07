import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import sys
import pdb

def readin_one_chromosome(QTL_dir, file_suffix):
    df = pd.DataFrame()
    for chromosome in range(1,2):
	fn1 = '%s/CHR%d_%s' % (QTL_dir, chromosome, file_suffix)
	try:
	    dfi = pd.read_csv(fn1, sep='\t')
	except:
	    print('Result for chr%d not exist' % chromosome)
	    continue
	df = df.append(dfi)
    df_sig = df.loc[df.groupby('PeakID')["P-value"].idxmin()]
    df_sig['BH_pvalue'] = multipletests(df_sig['P-value'], method='fdr_bh')[1]
    df_sig = df_sig[df_sig['BH_pvalue'] < 0.05]
    return [df[['PeakID', 'CHR_POS', 'P-value']], df_sig[['PeakID', 'CHR_POS', 'P-value', 'BH_pvalue']]]


def readin_both_results(QTL_dir, window, suffix):
    
    [real_QTL, real_QTL_sig] = readin_one_chromosome(QTL_dir, 'caQTLs_WINDOW_%skb%s.txt' % (str(window/1000.0), '_realGT_all'))
    [call_QTL, call_QTL_sig] = readin_one_chromosome(QTL_dir, 'caQTLs_WINDOW_%skb%s.txt' % (str(window/1000.0), suffix))
    [call_QTL_imputed, call_QTL_imputed_sig] = readin_one_chromosome('%s/Imputation' % QTL_dir, 'caQTLs_WINDOW_%skb_withImputation%s.txt' % (str(window/1000.0), suffix))

    call_QTL.columns = ['PeakID', 'CHR_POS', 'Pvalue_genotype_caller']
    call_QTL_sig.columns = ['PeakID', 'CHR_POS', 'Pvalue_genotype_caller', 'BH_pvalue_genotype_caller']
   
    call_QTL_imputed.columns = ['PeakID', 'CHR_POS', 'Pvalue_imputation']
    call_QTL_imputed_sig.columns = ['PeakID', 'CHR_POS', 'Pvalue_imputation', 'BH_pvalue_imputation']
 
    QTLs_all = real_QTL.merge(call_QTL, on = ['CHR_POS', 'PeakID'], how = 'outer').merge(call_QTL_imputed, on = ['CHR_POS', 'PeakID'], how = 'outer')
    QTLs_sig_from_all = real_QTL_sig.merge(call_QTL_sig, on = ['CHR_POS', 'PeakID'], how = 'outer').merge(call_QTL_imputed_sig, on = ['CHR_POS', 'PeakID'], how = 'outer')

    return [QTLs_all, QTLs_sig_from_all]



def obtain_metrics(peak_calling, window, suffix = '_noWeight'):

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
        
    [df, df_sig] = readin_both_results(QTL_dir, window, suffix = suffix)
    

    # genotype caller
    both_tested = df[(~df['Pvalue_genotype_caller'].isnull()) & (~df['P-value'].isnull())]
    r2_caller = np.corrcoef(-np.log10(both_tested['P-value']), -np.log10(both_tested['Pvalue_genotype_caller']))[0, 1]

    # imputation
    both_tested = df[(~df['Pvalue_imputation'].isnull()) & (~df['P-value'].isnull())]
    r2_imputation = np.corrcoef(-np.log10(both_tested['P-value']), -np.log10(both_tested['Pvalue_imputation']))[0, 1]

    print('Using genotype caller the Pearson correlation = %.2f; \n Using imputation the Pearson correlation = %.2f' % (r2_caller, r2_imputation))

    set1 = np.unique(df_sig[~df_sig['P-value'].isnull()].apply(lambda x: '_'.join((x['CHR_POS'], x['PeakID'])), axis=1))
    set2 = np.unique(df_sig[~df_sig['Pvalue_genotype_caller'].isnull()].apply(lambda x: '_'.join((x['CHR_POS'], x['PeakID'])), axis=1))
    set3 = np.unique(df_sig[~df_sig['Pvalue_imputation'].isnull()].apply(lambda x: '_'.join((x['CHR_POS'], x['PeakID'])), axis=1))

    precision_caller = len(np.intersect1d(set1, set2)) / float(len(set2))
    recall_caller = len(np.intersect1d(set1, set2)) / float(len(set1))
    precision_imputation = len(np.intersect1d(set1, set3) )/ float(len(set3))
    recall_imputation = len(np.intersect1d(set1, set3)) / float(len(set1))

    result = [peak_calling, window, r2_caller, r2_imputation, precision_caller, precision_imputation, recall_caller, recall_imputation]

    return result



if __name__ == '__main__':
    for WINDOW in [0, 1000]:
    	result = obtain_metrics('macs2', WINDOW)
	print(result)
