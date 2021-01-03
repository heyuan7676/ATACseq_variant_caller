import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import sys
sys.path.append('/work-zfs/abattle4/heyuan/old_work_files/yuan/tools/python_lib/lib/python2.7/site-packages')
import pdb


def readin_chromosomes_peakLevel(QTL_dir, file_suffix):

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




def readin_peakLevel_results(QTL_dir, window, suffix):

    real_QTL_sig = readin_chromosomes_peakLevel(QTL_dir, 'caQTLs_WINDOW_%skb%s.txt' % (str(window/1000.0), '_realGT_all'))
    call_QTL_sig = readin_chromosomes_peakLevel(QTL_dir, 'caQTLs_WINDOW_%skb%s.txt' % (str(window/1000.0), suffix))
    call_QTL_imputed_sig = readin_chromosomes_peakLevel('%s/Imputation' % QTL_dir, 'caQTLs_WINDOW_%skb_withImputation%s.txt' % (str(window/1000.0), suffix))

    call_QTL_sig.columns = ['PeakID', 'CHR_POS_called', 'Pvalue_genotype_caller', 'BH_pvalue_genotype_caller']
    call_QTL_imputed_sig.columns = ['PeakID', 'CHR_POS_imputation', 'Pvalue_imputation', 'BH_pvalue_imputation']

    QTLs_sig_from_all = real_QTL_sig.merge(call_QTL_sig, on = ['PeakID'], how = 'outer').merge(call_QTL_imputed_sig, on = ['PeakID'], how = 'outer')

    return QTLs_sig_from_all



# In[4]:


def readin_all_df(peak_calling, window, suffix = '_noWeight', readin_true_df = False):

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

    df_peakLevel = readin_peakLevel_results(QTL_dir, window, suffix = suffix)
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

    df_subset = append_R2_info(df_subset, col)
    return df_subset



def compare_QTLs(window):
    [df_pl, df_true_all] = readin_all_df(peak_calling, window, readin_true_df=True)
    N_called = np.sum(~df_pl['CHR_POS_called'].isnull())
    N_imputed = np.sum(~df_pl['CHR_POS_imputation'].isnull())
    N_true = np.sum(~df_pl['CHR_POS'].isnull())

    df_called = compare_top_SNPs(df_pl, df_true_all, col = 'CHR_POS_called')
    df_imputed = compare_top_SNPs(df_pl, df_true_all, col = 'CHR_POS_imputation')

    metrics = []
    for R2 in [0.2, 0.4, 0.6, 0.8, 1]:
        for THR in [1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 0]:
	    N_called = np.sum((df_called['df_called_pvalues_in_true'] - df_called['P-value'] <= THR) | df_called['R2'] >= R2)
	    N_imputed = np.sum((df_imputed['df_called_pvalues_in_true'] - df_imputed['P-value'] <= THR) | df_imputed['R2'] >= R2)
            precision_called = N_called / float(sum(~df_pl['Pvalue_genotype_caller'].isnull()))
            precision_imputed = N_imputed / float(sum(~df_pl['Pvalue_imputation'].isnull()))
    
            recall_called = N_called / float(sum(~df_pl['P-value'].isnull()))
            recall_imputed = N_imputed / float(sum(~df_pl['P-value'].isnull()))
    
            metrics.append([window, R2, THR, N_called, N_imputed, N_true, precision_called, precision_imputed, recall_called, recall_imputed])
    
    metrics = pd.DataFrame(metrics)
    metrics.columns = ['Window', 'R2', 'THR', 'N_called','N_imputed','N_true','precision_called', 'precision_imputed', 'recall_called', 'recall_imputed']
    print(metrics)
    return metrics
    

def append_R2_info(df_true_with_called, col):
    df_true_with_called_chromosome = np.array([int(x.split('_')[0]) for x in df_true_with_called[col]])
    df_true_with_called_withR2 = pd.DataFrame()

    for chromosome in range(1, 23):
        df_true_with_called_chr = df_true_with_called.iloc[np.where(df_true_with_called_chromosome == chromosome)[0]]
        df_true_with_called_chr = df_true_with_called_chr.merge(SNP_pair_r2[str(chromosome)], left_on=['CHR_POS', col], right_on=['CHR_A_BP_A', 'CHR_B_BP_B'], how = 'left')
    
        df_true_with_called_chr.loc[df_true_with_called_chr['CHR_POS'] == df_true_with_called_chr[col],'R2'] = 1
        df_true_with_called_chr.loc[df_true_with_called_chr['R2'].isnull(), 'R2'] = 0
        df_true_with_called_withR2 = df_true_with_called_withR2.append(df_true_with_called_chr)
        
    return df_true_with_called_withR2


if __name__ == '__main__':
    peak_calling = 'macs2'
    minDP = 2

    r2_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Genotype/plink'
    SNP_pair_r2 = dict()
    for chromosome in range(1, 23):
        SNP_pair_r2[str(chromosome)] = pd.read_csv('%s/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.r2.ld.variantPairs.txt' % (r2_dir, chromosome), delim_whitespace=True)


    metrics_peakLevel = pd.DataFrame()
    metrics_topVariant = pd.DataFrame()
    for w in [0, 1000, 10000, 100000]:
        [df, _] = readin_all_df('macs2', w)
        metrics_peakLevel = metrics_peakLevel.append(examine_peak_overlap(df, w))
	metrics_topVariant = metrics_topVariant.append(compare_QTLs(w))
    metrics_peakLevel.to_csv('Evaluation_metrics/%s_minDP%d_withImputation_noWeight_peakLevel_withR2.txt' % (peak_calling, minDP), index = False)
    metrics_topVariant.to_csv('Evaluation_metrics/%s_minDP%d_withImputation_noWeight_topVariant_withR2.txt' % (peak_calling, minDP), index = False)

