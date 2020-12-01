import pdb
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ranksums
from Evaluation_metrics import *
from prepare_genotype_data import *

def compute_genotype_metrics_called_and_imputed(sample, restrict_to_SNP = True):
    if restrict_to_SNP:
        suffix = '_allSNPs'
    else:
	suffix = '_allVariants'

    try:
	print('For %s:' % sample)
        # WGS genotype
        WGS_df = read_in_WGS_GT(sample)
    except:
	print('%s does not have genotype data' % s)	

    # variant calling information
    print('Evaluate genotype originally called from ATAC-seq reads...')
    orginally_called = obtain_atac_variants_df(sample, WGS_result=WGS_df, restrict_to_SNP=restrict_to_SNP, Imputed = False, return_metric = False)
    performance_original = compute_metric(orginally_called, sample = sample, minDP = minDP)

    print('Evalutae genotype imputed')
    imputed = obtain_atac_variants_df(sample, WGS_result=WGS_df, restrict_to_SNP=restrict_to_SNP, Imputed = True, return_metric = False)
    N_notExist_in_WGS = sum(imputed['REF_x'].isnull())

    print('Among variants presented in 1000 Genome project phase three')
    imputed = imputed[~imputed['REF_x'].isnull()]
    performance_imputed = compute_metric(imputed, sample = sample, minDP = 2)

    # compare the two called genotype datasets
    orginally_called = orginally_called[~orginally_called['REF_y'].isnull()] 
    imputed = imputed[~imputed['REF_y'].isnull()]

    combined = orginally_called.merge(imputed, on = ['#CHROM', 'POS'], how = 'outer')

    # categorize and evaluate
    only_in_orignal = combined[combined['REF_y_y'].isnull() & (~combined['REF_y_x'].isnull())]
    only_in_imputed = combined[combined['REF_y_x'].isnull() & (~combined['REF_y_y'].isnull())]
    called_in_both = combined[(~combined['REF_y_x'].isnull()) & (~combined['REF_y_y'].isnull())]
    assert len(called_in_both) + len(only_in_orignal) + len(only_in_imputed)== len(combined)

    # called by only one method
    only_in_orignal = only_in_orignal[only_in_orignal.columns[:7]] 
    only_in_orignal.columns = ['#CHROM', 'POS', 'REF_x', sample, 'REF_y', '%s_called' % sample, '%s_makeup' % sample]

    only_in_imputed = only_in_imputed[list(only_in_imputed.columns[:2]) + list(only_in_imputed.columns[7:])]
    only_in_imputed.columns = ['#CHROM', 'POS', 'REF_x', sample, 'REF_y', '%s_called' % sample, '%s_makeup' % sample]

    print('Evalutae genotype imputed - captured only in originally called variants')
    performance_only_in_original = compute_metric(only_in_orignal, sample = sample, minDP = 2)

    print('Evalutae genotype imputed - captured only in imputation')
    performance_only_in_imputed = compute_metric(only_in_imputed, sample = sample, minDP = 2) 

    performance_only_in_original.append('performance_only_in_original')
    performance_only_in_imputed.append('performance_only_in_imputed')

    performance = pd.DataFrame([performance_only_in_original, performance_only_in_imputed])
    performance.columns = ['Sample', 'minDP',
                           'N_variants_by_WGS', 'N_variants_by_ATAC', 'N_overlap_variants_by_ATAC', 'Recovered_percentage', 'Correctly_performance_percentage',
                           'Recall_AA', 'Recall_AB', 'Recall_BB', 'Precision_AA', 'Precision_AB', 'Precision_BB', 'group']

    # called by both
    print('Evalutae genotype imputed - captured in both methods')
    true_gt = map(convert_gt_to_number, np.array(called_in_both[['%s_x' % sample, '%s_called_x' % sample, '%s_called_y' % sample]]))
    true_gt = pd.DataFrame(np.array(true_gt))
    true_gt.columns = ['True', 'Originally_called', 'Imputed']
    true_gt['#CHROM'] = np.array(called_in_both['#CHROM'])
    true_gt['POS'] = np.array(called_in_both['POS'])

    discrepancy = true_gt[true_gt['Originally_called'] != true_gt['Imputed']]
    discrepancy = discrepancy.copy()
    discrepancy['CHR_POS'] = discrepancy[['#CHROM', 'POS']].apply(lambda x: '_'.join((str(x[0]), str(x[1]))), axis=1)

    stats = pd.DataFrame([sample, len(combined), N_notExist_in_WGS, len(only_in_orignal), len(called_in_both) - len(discrepancy), len(discrepancy), len(only_in_imputed)]).transpose()
    stats.columns = ["sample", "called_and_in_WGS", "called(imputed)_not_in_WGS","only_in_orignal", "called_in_both_consistent", "called_in_both_inconsistently", "only_in_imputed"]
    print(stats)

    print('Evalutae genotype imputed - captured in both methods but with disconcordant genotype')
    df = discrepancy.copy()
    original_precision_AB = sum((df['Originally_called'] == df['True']) & (df['Originally_called'] == 1)) / float(np.sum(df['Originally_called'] == 1))
    original_precision_AA_BB = sum((df['Originally_called'] == df['True']) & (df['Originally_called'] != 1)) / float(np.sum(df['Originally_called'] != 1))

    imputed_precision_AB = sum((df['Imputed'] == df['True']) & (df['Imputed'] == 1)) / float(np.sum(df['Imputed'] == 1))
    imputed_precision_AA_BB = sum((df['Imputed'] == df['True']) & (df['Imputed'] != 1)) / float(np.sum(df['Imputed'] != 1))

    discrepancy_metrics = pd.DataFrame([sample, original_precision_AB, original_precision_AA_BB, imputed_precision_AB, imputed_precision_AA_BB]).transpose()
    discrepancy_metrics.columns = ["sample", "original_precision_AB", "original_precision_AA_BB", "imputed_precision_AB", "imputed_precision_AA_BB"]
    print(discrepancy_metrics)

    if not os.path.exists('performance/combined'):
        os.makedirs('performance/combined')
    performance.to_csv('performance/combined/performance_one_method_%s%s.txt' % (sample, suffix), sep='\t', index = False)
    stats.to_csv('performance/combined/called_in_both_%s%s.txt' % (sample, suffix), sep='\t', index = False)
    discrepancy_metrics.to_csv('performance/combined/called_inconsistently_%s%s.txt' % (sample, suffix), sep='\t', index = False)


if __name__ == '__main__':
    minDP = 2
    sample = sys.argv[1]

    compute_genotype_metrics_called_and_imputed(sample, restrict_to_SNP = False)
