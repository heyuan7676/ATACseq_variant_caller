import numpy as np
import pandas as pd
import sys
from statsmodels.stats.multitest import multipletests
import pdb

def readin_fastQTL(method, peak_calling, cisDist, peakLevel = True):
    if peakLevel:
        df_peakLevel = pd.read_csv('%s/%s/%s/window%s.fastq.permutation.results.BH.txt' % (QTL_dir, method, peak_calling, cisDist), sep=' ')
        df_peakLevel = df_peakLevel[['pid', 'sid', 'npval', 'bh']]
        df_peakLevel.columns = ['peakID', 'SNP_%s' % method, 'QTL_level_Pvalue_%s' % method, 'BH_pvalue_%s' % method]
        return df_peakLevel
    else:
        df_QTLLevel = pd.DataFrame()
        for chromosome in range(1, 23):
            dfi = pd.read_csv('%s/%s/%s/chr%d.window%d.fastq.results' % (QTL_dir, method, peak_calling, chromosome, cisDist), sep=' ', header = None)
            df_QTLLevel = df_QTLLevel.append(dfi)
        df_QTLLevel.columns = ['peakID', 'SNP_%s' % method, 'distance', 'QTL_level_Pvalue_%s' % method, 'slope']
        return df_QTLLevel

   



def readin_called_QTLs(peak_calling, cisDist):

    QTL_real = readin_fastQTL('', peak_calling, cisDist, peakLevel = False) 

    QTL_gc = readin_fastQTL('VCF_files', peak_calling, cisDist)
    QTL_imputation = readin_fastQTL('Imputation', peak_calling, cisDist)
    QTL_integration = readin_fastQTL('Integration', peak_calling, cisDist)

    QTLs_all_peakLevel = QTL_real.merge(QTL_gc, left_on = ['peakID', 'SNP_'], right_on = ['peakID', 'SNP_VCF_files'], how = 'outer').merge(QTL_imputation, left_on = ['peakID', 'SNP_'], right_on =['peakID', 'SNP_Imputation'], how = 'outer').merge(QTL_integration, left_on = ['peakID', 'SNP_'], right_on = ['peakID', 'SNP_Integration'], how = 'outer')

    return QTLs_all_peakLevel


def readin_real_QTLs(peak_calling, cisDist):

    QTL_real = readin_fastQTL('', peak_calling, cisDist) 
    QTL_gc = readin_fastQTL('VCF_files', peak_calling, cisDist, peakLevel = False)
    QTL_imputation = readin_fastQTL('Imputation', peak_calling, cisDist, peakLevel = False)
    QTL_integration = readin_fastQTL('Integration', peak_calling, cisDist, peakLevel = False)

    QTLs_all_peakLevel = QTL_real.merge(QTL_gc, left_on = ['peakID', 'SNP_'], right_on = ['peakID', 'SNP_VCF_files'], how = 'outer').merge(QTL_imputation, left_on = ['peakID', 'SNP_'], right_on =['peakID', 'SNP_Imputation'], how = 'outer').merge(QTL_integration, left_on = ['peakID', 'SNP_'], right_on = ['peakID', 'SNP_Integration'], how = 'outer')

    return QTLs_all_peakLevel


def replicate_in_real(peak_calling):
    metrics_pi = []

    for w in window_list:
        df_pl = readin_called_QTLs(peak_calling, w)
        for method in ['VCF_files', 'Imputation', 'Integration']:
            col = 'BH_pvalue_%s' % method
            called_QTLs = df_pl[df_pl[col] < FDR]
            # set pvalue = 1 if not the top variant
            called_QTLs = called_QTLs.copy()
            called_QTLs.loc[called_QTLs['SNP_'].isnull(), 'QTL_level_Pvalue_'] = 1
            called_QTLs['BH_replication'] = multipletests(called_QTLs['QTL_level_Pvalue_'], method = 'fdr_bh')[1]
            true_positive = np.sum(called_QTLs['BH_replication'] < FDR)
            pi1 = true_positive / float(len(called_QTLs))
            print([w, method, pi1])
            metrics_pi.append([peak_calling, w, method, pi1])

    metrics_pi = pd.DataFrame(metrics_pi)
    metrics_pi.columns = ['peak_calling', 'window', 'method', 'pi1']
    metrics_pi.to_csv('Evaluation_metrics/%s_minDP%d_topVariant_cisDist_pi1_called_repInReal.txt' % (peak_calling.replace('/', '_'), minDP), sep='\t', index = False)
 

def replicate_in_called(peak_calling):
    metrics_pi = []

    for w in window_list:
        df_pl = readin_real_QTLs(peak_calling, w)
        called_QTLs = df_pl[df_pl['BH_pvalue_'] < FDR]
        called_QTLs = called_QTLs.copy()
    
        for method in ['VCF_files', 'Imputation', 'Integration']:
            called_QTLs.loc[called_QTLs['SNP_%s' % method].isnull(), 'QTL_level_Pvalue_%s' % method] = 1
            called_QTLs['BH_replication'] = multipletests(called_QTLs['QTL_level_Pvalue_%s' % method], method = 'fdr_bh')[1]
            true_positive = np.sum(called_QTLs['BH_replication'] < FDR)
            pi1 = true_positive / float(len(called_QTLs))
            print([w, method, pi1])
            metrics_pi.append([peak_calling, w, method, pi1])

    metrics_pi = pd.DataFrame(metrics_pi)
    metrics_pi.columns = ['peak_calling', 'window', 'method', 'pi1']
    metrics_pi.to_csv('Evaluation_metrics/%s_minDP%d_topVariant_cisDist_pi1_real_repInCalled.txt' % (peak_calling.replace('/', '_'), minDP), sep='\t', index = False)




if __name__ == '__main__':
    peak_calling = sys.argv[1]
    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie'
    minDP = 3
    window_list = [300, 500, 1000, 10000, 100000, 1000000]

    FDR = 0.05
    QTL_dir = '%s/fastQTL' % (root_dir)

    replicate_in_real(peak_calling)
    replicate_in_called(peak_calling)
