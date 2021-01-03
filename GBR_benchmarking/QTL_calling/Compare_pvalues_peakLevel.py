import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import sys
sys.path.append('/work-zfs/abattle4/heyuan/old_work_files/yuan/tools/python_lib/lib/python2.7/site-packages')
import pdb


def readin_chromosomes(QTL_dir, file_suffix):
    df = pd.DataFrame()
    for chromosome in range(1,23):
        fn1 = '%s/CHR%d_%s' % (QTL_dir, chromosome, file_suffix)
        try:
            dfi = pd.read_csv(fn1, sep='\t')
        except:
            print('Result for chr%d not exist' % chromosome)
            continue
        df = df.append(dfi)

    df = df[~df['P-value'].isnull()]
    df = df[~df[['PeakID', 'CHR_POS']].duplicated()].reset_index(drop = True)
    
    return df[['PeakID', 'CHR_POS', 'P-value', 'Start', 'End', "POS"]]




def readin_peakLevel_results(QTL_dir, window):

    real_QTL_sig = compute_QTL_peakLevel(window, QTL_dir, saveSuffix = '_realGT_all')
    call_QTL_sig = compute_QTL_peakLevel(window, QTL_dir, saveSuffix = ''))
    call_QTL_imputed_sig = compute_QTL_peakLevel(window, '%s/Imputation' % QTL_dir, saveSuffix = '_Imputation'))

    call_QTL_sig.columns = ['PeakID', 'CHR_POS_called', 'Pvalue_genotype_caller', 'BH_pvalue_genotype_caller']
    call_QTL_imputed_sig.columns = ['PeakID', 'CHR_POS_imputation', 'Pvalue_imputation', 'BH_pvalue_imputation']

    QTLs_sig_from_all = real_QTL_sig.merge(call_QTL_sig, on = ['PeakID'], how = 'outer').merge(call_QTL_imputed_sig, on = ['PeakID'], how = 'outer')

    return QTLs_sig_from_all



# In[4]:


def readin_all_df(peak_calling, window, readin_peak_level = True, readin_true_df = False):

    alignment_dir = 'alignment_bowtie' # alignment_subsample_0.5
    GT_subDir = 'minDP%d' % minDP

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

    if readin_peak_level:
        df_peakLevel = readin_peakLevel_results(QTL_dir, window)
    else:
	df_peakLevel = None

    if readin_true_df:
        true_df_all = readin_chromosomes(QTL_dir, 'caQTLs_WINDOW_%skb%s.txt' % (str(window/1000.0), '_realGT_all'))
    else:
        true_df_all = None
    
    return [df_peakLevel, true_df_all]


# In[5]:


# Peak Overlap
def examine_peak_overlap(df_peakLevel, window):
    set1 = np.unique(df_peakLevel[~df_peakLevel['P-value'].isnull()]['PeakID'])
    set2 = np.unique(df_peakLevel[~df_peakLevel['Pvalue_genotype_caller'].isnull()]['PeakID'])
    set3 = np.unique(df_peakLevel[~df_peakLevel['Pvalue_imputation'].isnull()]['PeakID'])

    precision_caller = len(np.intersect1d(set1, set2)) / float(len(set2))
    recall_caller = len(np.intersect1d(set1, set2)) / float(len(set1))

    precision_imputation = len(np.intersect1d(set1, set3) )/ float(len(set3))
    recall_imputation = len(np.intersect1d(set1, set3)) / float(len(set1))

    result = [window, len(set1), len(set2), len(set3), precision_caller, precision_imputation, recall_caller, recall_imputation]
    result = pd.DataFrame(result).transpose()
    result.columns = ["window", "Number_QTL_true", "Number_QTL_gc", "Number_QTL_imputation", "precision_caller", "precision_imputation", "recall_caller", "recall_imputation"]
   
    print(result)
 
    return(result)




def compare_top_SNPs(df_peakLevel_sig, df_true_all, col = 'CHR_POS_called'):
    df_subset = df_peakLevel_sig[['PeakID','CHR_POS', col, 'P-value']]
    df_subset = df_subset[~df_subset[col].isnull()]
  
    df_true_pvalue = df_subset.merge(df_true_all, left_on = ['PeakID', col], right_on = ['PeakID', 'CHR_POS'], how = 'left') 
    df_subset = df_subset.copy()
    df_subset['df_called_pvalues_in_true'] = np.array(df_true_pvalue['P-value_y'])
    df_subset.loc[df_subset['df_called_pvalues_in_true'].isnull(), 'df_called_pvalues_in_true'] = 100
    df_subset.loc[df_subset['P-value'].isnull(), 'P-value'] = -100
    return df_subset



def compare_QTLs(window, df_true_all_for_recall, N_true_for_recall):
    [df_pl, df_true_all] = readin_all_df(peak_calling, window, readin_true_df=True)
    N_called = np.sum(~df_pl['CHR_POS_called'].isnull())
    N_imputed = np.sum(~df_pl['CHR_POS_imputation'].isnull())
    N_true = np.sum(~df_pl['CHR_POS'].isnull())

    df_called = compare_top_SNPs(df_pl, df_true_all, col = 'CHR_POS_called')
    df_imputed = compare_top_SNPs(df_pl, df_true_all, col = 'CHR_POS_imputation')

    metrics = []
    for THR in THR_list:
        precision_called = np.sum(df_called['df_called_pvalues_in_true'] - df_called['P-value'] <= THR) / float(sum(~df_pl['Pvalue_genotype_caller'].isnull()))
        precision_imputed = np.sum(df_imputed['df_called_pvalues_in_true'] - df_imputed['P-value'] <= THR) / float(sum(~df_pl['Pvalue_imputation'].isnull()))
    
        metrics.append([window, THR, N_called, N_imputed, N_true, precision_called, precision_imputed, 'Precision'])
   
    # recall
    df_called = compare_top_SNPs(df_pl, df_true_all_for_recall, col = 'CHR_POS_called')
    df_imputed = compare_top_SNPs(df_pl, df_true_all_for_recall, col = 'CHR_POS_imputation') 
    for THR in THR_list:
        recall_called = np.sum((df_called['df_called_pvalues_in_true'] - df_called['P-value']) <= THR) / float(N_true_for_recall)
        recall_imputed = np.sum((df_imputed['df_called_pvalues_in_true'] - df_imputed['P-value']) <= THR) / float(N_true_for_recall)

	metrics.append([window, THR, N_called, N_imputed, N_true, precision_called, precision_imputed, 'Recall'])
 
    metrics = pd.DataFrame(metrics)
    metrics.columns = ['Window', 'THR', 'N_called','N_imputed','N_true','Called', 'Imputation', 'Metric']

    print(metrics)
    return metrics
    


if __name__ == '__main__':
    peak_calling = 'macs2'
    minDP = 2
    metrics_peakLevel = pd.DataFrame()
    metrics_topVariant = pd.DataFrame()

    THR_list = [1e-3, 1e-5, 1e-10, 0]
    w = 10000
    [df, all_qtl_tests] = readin_all_df('macs2', w, readin_true_df=True)
    N_all = np.sum(~df['CHR_POS'].isnull())
    metrics_peakLevel = metrics_peakLevel.append(examine_peak_overlap(df, w))
    metrics_topVariant = metrics_topVariant.append(compare_QTLs(w, all_qtl_tests, N_all))

    for w in [0, 1000, 10000, 10000][:2]:
        [df, _] = readin_all_df('macs2', w)
        metrics_peakLevel = metrics_peakLevel.append(examine_peak_overlap(df, w))
	metrics_topVariant = metrics_topVariant.append(compare_QTLs(w, all_qtl_tests, N_all))
    metrics_peakLevel.to_csv('Evaluation_metrics/%s_minDP%d_withImputation_noWeight_peakLevel.txt' % (peak_calling, minDP), index = False)
    metrics_topVariant.to_csv('Evaluation_metrics/%s_minDP%d_withImputation_noWeight_topVariant.txt' % (peak_calling, minDP), index = False)
