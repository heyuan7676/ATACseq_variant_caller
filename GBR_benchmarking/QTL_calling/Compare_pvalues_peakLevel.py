import numpy as np
import pandas as pd
import sys
from statsmodels.stats.multitest import multipletests
import pdb

def readin_chromosomes(QTL_dir, minDP, cisDist, method):
    try:
        df_peakLevel = pd.read_csv('%s/matrixeQTL_minDP%d_cisDist_%dkb_%s_peakLevel.txt' % (QTL_dir, minDP, (cisDist)/1000, method), sep='\t')
        df_peakLevel.columns = ['gene', 'SNP_%s' % method, 'QTL_level_Pvalue_%s' % method, 'Peak_level_Pvalue_%s' % method, 'BH_pvalue_%s' % method]
        return df_peakLevel
    
    except:
        print('Compute peal-level statistics')
        
    df_peakLevel = pd.DataFrame()
    for chromosome in range(1,23):
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
            
        df_peakLevel = df_peakLevel.append(dfi_peakLevel)
        
    df_peakLevel['BH_pvalue'] = multipletests(df_peakLevel['p-value'], method='fdr_bh')[1]
    df_peakLevel = df_peakLevel[['gene', 'SNP', 'p-value', 'Peak_level_Pvalue', 'BH_pvalue']]
    df_peakLevel.columns = ['gene', 'SNP_%s' % method, 'QTL_level_Pvalue_%s' % method,'Peak_level_Pvalue_%s' % method, 'BH_pvalue_%s' % method]
    
    df_peakLevel.to_csv('%s/matrixeQTL_minDP%d_cisDist_%dkb_%s_peakLevel.txt' % (QTL_dir, minDP, (cisDist)/1000, method), sep='\t', index=False)

    return df_peakLevel


def readin_all_results(QTL_dir,  minDP, cisDist):

    QTL_real = readin_chromosomes(QTL_dir, minDP, cisDist, 'WGS')
    QTL_gc = readin_chromosomes(QTL_dir, minDP, cisDist, 'GC')
    QTL_imputation = readin_chromosomes(QTL_dir, minDP, cisDist, 'Imputation')
    QTL_integration = readin_chromosomes(QTL_dir, minDP, cisDist, 'Integration')

    QTLs_all_peakLevel = QTL_real.merge(QTL_gc, on = ['gene'], how = 'outer').merge(QTL_imputation, on = ['gene'], how = 'outer').merge(QTL_integration, on = ['gene'], how = 'outer')

    return QTLs_all_peakLevel



def append_R2_info(df_true_with_called):
    df_true_with_called_chromosome = np.array([int(x.split('_')[0]) for x in df_true_with_called['SNP_compare']])
    df_true_with_called_withR2 = pd.DataFrame()

    for chromosome in range(1, 23):
        df_true_with_called_chr = df_true_with_called.iloc[np.where(df_true_with_called_chromosome == chromosome)[0]]
        df_true_with_called_chr = df_true_with_called_chr.merge(SNP_pair_r2[str(chromosome)], left_on=['SNP_WGS', 'SNP_compare'], right_on=['CHR_A_BP_A', 'CHR_B_BP_B'], how = 'left')
        df_true_with_called_chr.loc[df_true_with_called_chr['R2'].isnull(), 'R2'] = 0
        df_true_with_called_withR2 = df_true_with_called_withR2.append(df_true_with_called_chr)

    return df_true_with_called_withR2



# Peak Overlap
def examine_peak_overlap(df_peakLevel, window):
    set1 = np.unique(df_peakLevel[~df_peakLevel['BH_pvalue_WGS'].isnull()]['gene'])
    set2 = np.unique(df_peakLevel[~df_peakLevel['BH_pvalue_GC'].isnull()]['gene'])
    set3 = np.unique(df_peakLevel[~df_peakLevel['BH_pvalue_Imputation'].isnull()]['gene'])
    set4 = np.unique(df_peakLevel[~df_peakLevel['BH_pvalue_Integration'].isnull()]['gene'])

    precision_caller = len(np.intersect1d(set1, set2)) / float(len(set2))
    recall_caller = len(np.intersect1d(set1, set2)) / float(len(set1))

    precision_imputation = len(np.intersect1d(set1, set3) )/ float(len(set3))
    recall_imputation = len(np.intersect1d(set1, set3)) / float(len(set1))
    
    precision_integration = len(np.intersect1d(set1, set4) )/ float(len(set4))
    recall_integration = len(np.intersect1d(set1, set4)) / float(len(set1))

    #venn3([set(set1), set(set2), set(set3)], ('From Genotype data', 'Use variants from variant caller', 'Add variants from imputation'))
    #plt.show()
    #plt.close()

    result = [window, len(set1), len(set2), len(set3), len(set4), precision_caller, precision_imputation,precision_integration, recall_caller, recall_imputation, recall_integration]
    result = pd.DataFrame(result).transpose()
    result.columns = ["window", "Number_QTL_true", "Number_QTL_gc", "Number_QTL_imputation", "Number_QTL_integration","precision_caller", "precision_imputation", "precision_integration","recall_caller", "recall_imputation", "recall_integration"]
    
    return(result)




def readin_all_df(cisDist, method = 'WGS'):
    try:
        df = pd.read_csv('%s/matrixeQTL_minDP%d_cisDist_%dkb_%s.txt' % (QTL_dir, minDP, (cisDist)/1000, method), sep='\t')
        return df
    except:
        print('COmbine results from chromosomes')
        
    df = pd.DataFrame()
    for chromosome in range(1,23):
        fn1 = '%s/matrixeQTL_chromosome%d_minDP%d_cisDist_%dkb_%s.txt' % (QTL_dir, chromosome, minDP, (cisDist)/1000, method)
        try:
            dfi = pd.read_csv(fn1, sep='\t')
        except:
            print(fn1)
            print('Result for chr%d not exist' % chromosome)
            continue
            
        df = df.append(dfi)
        
    df['BH_pvalue'] = multipletests(df['p-value'], method='fdr_bh')[1]
    df = df[['gene', 'SNP', 'p-value', 'BH_pvalue']]
    df.columns = ['gene', 'SNP', 'QTL_level_Pvalue_%s' % method, 'BH_pvalue_%s' % method]
    
    df.to_csv('%s/matrixeQTL_minDP%d_cisDist_%dkb_%s.txt' % (QTL_dir, minDP, (cisDist)/1000, method), sep='\t', index=False)

    return df


# In[6]:


def compare_top_SNPs(df_peakLevel, df_true_all, col = ['SNP_called', 'Pvalue_genotype_caller']):

    df_subset = df_peakLevel[['gene','SNP_WGS', 'QTL_level_Pvalue_WGS']  + col]
    df_subset.columns = ['gene', 'SNP_WGS', 'Pvalue_WGS','SNP_compare', 'P-value_compare']
    df_subset = df_subset[~df_subset['P-value_compare'].isnull()]

    df_true_pvalue = df_subset.merge(df_true_all, left_on = ['gene', 'SNP_compare'], right_on = ['gene', 'SNP'], how = 'left')

    df_subset = df_subset.copy()
    df_subset['df_called_pvalues_in_true'] = np.array(df_true_pvalue['QTL_level_Pvalue_WGS'])
    df_subset.loc[df_subset['df_called_pvalues_in_true'].isnull(), 'df_called_pvalues_in_true'] = 100
    df_subset.loc[df_subset['Pvalue_WGS'].isnull(), 'Pvalue_WGS'] = -100

    return df_subset




def compare_QTLs(window, df_pl, df_true_all_for_recall, N_true_for_recall):
    
    N_true = np.sum(df_pl['BH_pvalue_WGS'] < FDR)
    metrics = []

    for method in ['GC', 'Imputation', 'Integration']:
        df_pl_method = df_pl[df_pl['BH_pvalue_%s' % method] < FDR]
        N_called = np.sum(df_pl['BH_pvalue_%s' % method] < FDR)

        df_called = compare_top_SNPs(df_pl_method, df_true_all_for_recall, col = ['SNP_%s' % method, 'QTL_level_Pvalue_%s' % method])
        df_called_withR2 = append_R2_info(df_called)

        for THR in THR_list:
            precision_called = np.sum(df_called['df_called_pvalues_in_true'] - df_called['Pvalue_WGS'] <= THR) / float(N_called)
            metrics.append([window, THR, method, N_called, N_true, precision_called, 'Precision', 'Pvalue_threshold'])

        for R2 in R2_list:
            precision_called = np.sum(df_called_withR2['R2']  >= R2) / float(N_called)
            metrics.append([window, R2, method, N_called, N_true, precision_called, 'Precision', 'R2_threshold'])

        # recall
        for THR in THR_list:
            recall_called = np.sum((df_called['df_called_pvalues_in_true'] - df_called['Pvalue_WGS']) <= THR) / float(N_true_for_recall)
            metrics.append([window, THR, method, N_called, N_true, recall_called,'Recall', 'Pvalue_threshold'])

        for R2 in R2_list:
            recall_called = np.sum(df_called_withR2['R2'] >= R2) / float(N_true_for_recall)
            metrics.append([window, R2, method, N_called, N_true, recall_called,'Recall', 'R2_threshold'])


    metrics = pd.DataFrame(metrics)
    metrics.columns = ['Window', 'THR', 'method','N_called','N_true','Value','Metric', 'Standard']

    print(metrics)
    return metrics



if __name__ == '__main__':
    peak_calling = sys.argv[1]
    minDP = int(sys.argv[2])
    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie'
    GT_subDir = ''

    FDR = 0.05

    if peak_calling == 'macs2':
        PEAK_dir = '%s/Peaks_MACS2' % root_dir
        QTL_dir = '%s/QTLs_MACS2/%s' % (root_dir, GT_subDir)
        WGS_QTL_dir = '%s/QTLs_MACS2/%s' % (root_dir, 'WGS')

    if peak_calling == 'macs2_combined':
        PEAK_dir = '%s/Peak_version1/Peaks/combined' % root_dir
        QTL_dir = '%s/QTLs_MACS2_combined/%s' % (root_dir, GT_subDir)
        WGS_QTL_dir = '%s/QTLs_MACS2_combined/%s' % (root_dir, 'WGS')

    elif peak_calling == 'Genrich':
        PEAK_dir = '%s/Peaks_Genrich' % root_dir
        QTL_dir = '%s/QTLs_Genrich/%s' % (root_dir, GT_subDir)
        WGS_QTL_dir = '%s/QTLs_Genrich/%s' % (root_dir, 'WGS')

    elif peak_calling == 'Genrich_combined':
        PEAK_dir = '%s/Peaks_Genrich/combined' % root_dir
        QTL_dir = '%s/QTLs_Genrich_combined/%s' % (root_dir, GT_subDir)
        WGS_QTL_dir = '%s/QTLs_Genrich_combined/%s' % (root_dir, 'WGS')


    r2_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/plink/'
    SNP_pair_r2 = dict()

    for chromosome in range(1, 23):
        print('chromosome%d' % chromosome)
        SNP_pair_r2[str(chromosome)] = pd.read_csv('%s/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.r2.ld.variantPairs.txt' % (r2_dir, chromosome), delim_whitespace=True)

    metrics_peakLevel = pd.DataFrame()
    metrics_topVariant = pd.DataFrame()

    QTL_all = readin_chromosomes(QTL_dir, minDP, 0, 'WGS')
    N_all = np.sum(QTL_all['BH_pvalue_WGS'] < FDR)

    THR_list = [1e-3, 1e-5, 1e-10, 0]
    R2_list = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7 ,0.8, 0.9, 1]
    for w in [0, 1000, 10000]:
        all_qtl_tests = readin_all_df(w)
        df_pl = readin_all_results(QTL_dir,  minDP, w)
        all_qtl_tests = all_qtl_tests.set_index('gene').loc[df_pl[df_pl['BH_pvalue_WGS'] < FDR]['gene']].reset_index()

        metrics_peakLevel = metrics_peakLevel.append(examine_peak_overlap(df_pl, w))
        metrics_topVariant = metrics_topVariant.append(compare_QTLs(w, df_pl, all_qtl_tests, N_all))

    metrics_peakLevel.to_csv('Evaluation_metrics/%s_minDP%d_peakLevel.txt' % (peak_calling, minDP), index = False)
    metrics_topVariant.to_csv('Evaluation_metrics/%s_minDP%d_topVariant.txt' % (peak_calling, minDP), index = False)


