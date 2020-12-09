import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import sys
import pdb

def readin_chromosomes(QTL_dir, file_suffix):

    try:
	df_peakLevel = pd.read_csv('%s/All.CHR_%s' % (QTL_dir, file_suffix), sep='\t')
	df_peakLevel = df_peakLevel[df_peakLevel['BH_pvalue'] < 0.05]
        return df_peakLevel[['PeakID', 'CHR_POS', 'P-value', 'BH_pvalue']]

    except:
	print('Read in QTL data')

    df_peakLevel = pd.DataFrame()
    for chromosome in range(1,23):
	fn1 = '%s/CHR%d_%s' % (QTL_dir, chromosome, file_suffix)
	try:
	    dfi_peakLevel = pd.read_csv('%s_peakLevel' % fn1, sep='\t')
	except:
	    dfi = pd.read_csv(fn1, sep='\t')
	    dfi = dfi[~dfi['P-value'].isnull()].reset_index(drop = True)
    	    dfi_peakLevel = dfi.loc[dfi.groupby('PeakID')["P-value"].idxmin()]
	    dfi_peakLevel.to_csv('%s_peakLevel' % fn1, sep='\t', index = False)
	df_peakLevel = df_peakLevel.append(dfi_peakLevel)

    df_peakLevel = df_peakLevel[['PeakID','CHR_POS', 'P-value', 'Start','End']]
    df_peakLevel['BH_pvalue'] = multipletests(df_peakLevel['P-value'], method='fdr_bh')[1]

    df_peakLevel = df_peakLevel.sort_values('BH_pvalue')
    df_peakLevel.to_csv('%s/All.CHR_%s' % (QTL_dir, file_suffix), sep='\t', index = False)

    df_peakLevel = df_peakLevel[df_peakLevel['BH_pvalue'] < 0.05]
    return df_peakLevel[['PeakID', 'CHR_POS', 'P-value', 'BH_pvalue']]


def readin_both_results(QTL_dir, window, suffix):
   
    real_QTL_sig = readin_chromosomes(QTL_dir, 'caQTLs_WINDOW_%skb%s.txt' % (str(window/1000.0), '_realGT_all'))
    call_QTL_sig = readin_chromosomes(QTL_dir, 'caQTLs_WINDOW_%skb%s.txt' % (str(window/1000.0), suffix))
    call_QTL_imputed_sig = readin_chromosomes('%s/Imputation' % QTL_dir, 'caQTLs_WINDOW_%skb_withImputation%s.txt' % (str(window/1000.0), suffix))

    call_QTL_sig.columns = ['PeakID', 'CHR_POS', 'Pvalue_genotype_caller', 'BH_pvalue_genotype_caller']
    call_QTL_imputed_sig.columns = ['PeakID', 'CHR_POS', 'Pvalue_imputation', 'BH_pvalue_imputation']
 
    QTLs_sig_from_all = real_QTL_sig.merge(call_QTL_sig, on = ['CHR_POS', 'PeakID'], how = 'outer').merge(call_QTL_imputed_sig, on = ['CHR_POS', 'PeakID'], how = 'outer')

    return QTLs_sig_from_all



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
        
    df_peakLevel = readin_both_results(QTL_dir, window, suffix = suffix)
    

    # genotype caller
    both_tested = df_peakLevel[(~df_peakLevel['Pvalue_genotype_caller'].isnull()) & (~df_peakLevel['P-value'].isnull())]
    r2_caller = np.corrcoef(-np.log10(both_tested['P-value']), -np.log10(both_tested['Pvalue_genotype_caller']))[0, 1]

    # imputation
    both_tested = df_peakLevel[(~df_peakLevel['Pvalue_imputation'].isnull()) & (~df_peakLevel['P-value'].isnull())]
    r2_imputation = np.corrcoef(-np.log10(both_tested['P-value']), -np.log10(both_tested['Pvalue_imputation']))[0, 1]

    print('Using genotype caller the Pearson correlation = %.2f; \n Using imputation the Pearson correlation = %.2f' % (r2_caller, r2_imputation))

    set1 = np.unique(df_peakLevel[~df_peakLevel['P-value'].isnull()].apply(lambda x: '_'.join((x['CHR_POS'], x['PeakID'])), axis=1))
    set2 = np.unique(df_peakLevel[~df_peakLevel['Pvalue_genotype_caller'].isnull()].apply(lambda x: '_'.join((x['CHR_POS'], x['PeakID'])), axis=1))
    set3 = np.unique(df_peakLevel[~df_peakLevel['Pvalue_imputation'].isnull()].apply(lambda x: '_'.join((x['CHR_POS'], x['PeakID'])), axis=1))

    precision_caller = len(np.intersect1d(set1, set2)) / float(len(set2))
    recall_caller = len(np.intersect1d(set1, set2)) / float(len(set1))

    precision_imputation = len(np.intersect1d(set1, set3) )/ float(len(set3))
    recall_imputation = len(np.intersect1d(set1, set3)) / float(len(set1))

    result = [peak_calling, window, r2_caller, r2_imputation, len(set1), len(set2), len(set3), precision_caller, precision_imputation, recall_caller, recall_imputation]

    return result



if __name__ == '__main__':
    WINDOW = int(sys.argv[1])
    peak_calling = 'macs2'
    result = obtain_metrics(peak_calling, WINDOW)
    result = pd.DataFrame(result).transpose()
    result.columns = ["peak_calling", "window", "r2_caller", "r2_imputation", "Number_QTL_true", "Number_QTL_gc", "Number_QTL_imputation", "precision_caller", "precision_imputation", "recall_caller", "recall_imputation"]
    print(result)
    #result.to_csv('Evaluation_metrics/%s_minDP2_%skb_withImputation_noWeight.txt' % (peak_calling, str(WINDOW/1000.0)), sep='\t')



