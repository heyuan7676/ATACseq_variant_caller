import pdb
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ranksums
from Evaluation_metrics import *


def convert_gt_to_number(arri):
    '''
    convert genotype type A/A to 0, 1, 2 ...
    '''
    arri_idx = np.where(arri!='0')[0]
    nts = np.unique([a for b in [x.split('/') for x in arri[arri_idx]] for a in b])
    # remove the sites with more than two alleles
    if len(nts) > 2:
        return np.ones(len(arri)) * (-1)
    nts_dict = {}
    i = 0
    for ni in nts:
        nts_dict[ni] = i
        i += 1
    numeric_gt = [np.sum([nts_dict[a] for a in x.split('/')]) for x in arri[arri_idx]]
    numeric_gt_arr = np.ones(len(arri)) * (-1)
    numeric_gt_arr[arri_idx] = numeric_gt
    return numeric_gt_arr



def compute_sens_recall(dat, sample, denominator_arr, group_text):
    dat.columns = ['#CHROM', 'POS', 'REF_x', sample, 'REF_y', '%s_called' % sample, '%s_makeup' % sample]
    df = obtain_confusion_matrix(dat, sample = sample)
    metric_list = [[sample]]
    metric_list.append(list(np.diag(df)))
    metric_list.append(list(np.diag(df) / np.sum(df, axis=0)))
    metric_list.append(list(np.diag(df) / np.sum(df, axis=1)))
    metric_list.append(list(denominator_arr))
    metric_list.append(list(np.diag(df) / denominator_arr))
    metric_list.append([group_text])

    metric_list = [a for b in metric_list for a in b]

    return metric_list 




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


    # read in variants from 1k Genome
    oneK_snps = pd.DataFrame()
    for c in range(1, 23):
	snps = read_in_1k_variants(c)
	oneK_snps = oneK_snps.append(snps)

    # variant calling information
    print('Evaluate genotype originally called from ATAC-seq reads...')
    orginally_called = obtain_atac_variants_df(sample, oneK_snps, WGS_result=WGS_df, restrict_to_SNP=restrict_to_SNP, Imputed = False)

    print('Evalutae genotype imputed')
    imputed = obtain_atac_variants_df(sample, oneK_snps, WGS_result=WGS_df, restrict_to_SNP=restrict_to_SNP, Imputed = True)
    N_notExist_in_WGS = sum(imputed['REF_x'].isnull())

    # compare the two called genotype datasets

    print('Among variants presented in 1000 Genome project phase three')
    #orginally_called = orginally_called[~orginally_called['REF_y'].isnull()] 
    #imputed = imputed[~imputed['REF_y'].isnull()]

    combined = orginally_called.merge(imputed, on = ['#CHROM', 'POS'], how = 'outer')
    combined = combined[~combined['REF_x_x'].isnull()]

    # get number of AA, AB, BB from the genotype data
    total_number_df = combined[list(combined.columns[:4]) + list(combined.columns[2:4]) + list(combined.columns[6:7])]
    total_number_df.columns = ['#CHROM', 'POS', 'REF_x', sample, 'REF_y', '%s_called' % sample, '%s_makeup' % sample]
    total_number_mx = obtain_confusion_matrix(total_number_df, sample)
    total_number_array = map(float, np.diag(total_number_mx)) 


    # categorize and evaluate
    in_orignal = combined[(~combined['REF_y_x'].isnull())]
    in_orignal = in_orignal[in_orignal.columns[:7]]

    in_imputed = combined[(~combined['REF_y_y'].isnull())]
    in_imputed = in_imputed[list(in_imputed.columns[:2]) + list(in_imputed.columns[7:])]

    print('Evalutae genotype imputed - captured in originally called variants')
    performance_in_original = compute_sens_recall(in_orignal, sample, total_number_array, 'performance_in_original')

    print('Evalutae genotype imputed - captured in imputation')
    performance_in_imputed = compute_sens_recall(in_imputed, sample, total_number_array, 'performance_in_imputed')


    # called by only one method

    print('Evalutae genotype imputed - captured only in originally called variants')
    only_in_orignal = combined[combined['REF_y_y'].isnull() & (~combined['REF_y_x'].isnull())]
    only_in_orignal = only_in_orignal[only_in_orignal.columns[:7]]
    performance_only_in_original = compute_sens_recall(only_in_orignal, sample, total_number_array, 'performance_only_in_original')

    print('Evalutae genotype imputed - captured only in imputation')
    only_in_imputed = combined[combined['REF_y_x'].isnull() & (~combined['REF_y_y'].isnull())]
    only_in_imputed = only_in_imputed[list(only_in_imputed.columns[:2]) + list(only_in_imputed.columns[7:])]
    performance_only_in_imputed = compute_sens_recall(only_in_imputed, sample, total_number_array, 'performance_only_in_imputed') 

    # called by both
    print('Evalutae genotype imputed - captured in both methods')
    called_in_both = combined[(~combined['REF_y_x'].isnull()) & (~combined['REF_y_y'].isnull())]
    called_in_both = called_in_both.copy()
    called_in_both['CHR_POS'] = called_in_both[['#CHROM', 'POS']].apply(lambda x: '_'.join((str(x[0]), str(x[1]))), axis=1)

    true_gt = map(convert_gt_to_number, np.array(called_in_both[['%s_x' % sample, '%s_called_x' % sample, '%s_called_y' % sample]]))
    true_gt = np.array(true_gt)
    called_in_both['True'] = true_gt[:, 0]
    called_in_both['Originally_called'] = true_gt[:, 1]
    called_in_both['Imputed'] = true_gt[:, 2]

    print('Evalutae genotype imputed - captured consistently from two sources')
    con_df = called_in_both[called_in_both['Originally_called'] == called_in_both['Imputed']]
    con_df = con_df[con_df.columns[:7]]
    performance_consistent = compute_sens_recall(con_df, sample, total_number_array, 'performance_consistent')

    print('Evalutae genotype imputed - captured In-consistently from two sources')
    discrepancy = called_in_both[called_in_both['Originally_called'] != called_in_both['Imputed']]
    discrepancy_use_original = discrepancy[discrepancy.columns[:7]]
    performance_discrepancy_use_original = compute_sens_recall(discrepancy_use_original, sample, total_number_array, 'performance_discrepancy_use_original')

    discrepancy_use_imputed = discrepancy[list(discrepancy.columns[:2]) + list(discrepancy.columns[7:12])]
    performance_discrepancy_use_imputed = compute_sens_recall(discrepancy_use_imputed, sample, total_number_array, 'performance_discrepancy_use_imputed')


    # save the performance
    performance = pd.DataFrame([performance_in_original, performance_in_imputed, performance_only_in_original, performance_only_in_imputed, performance_consistent, performance_discrepancy_use_original, performance_discrepancy_use_imputed])
    performance.columns = ['Sample', 'Number_AA_called_true', 'Number_AB_called_true', 'Number_BB_called_true', 'Precision_AA', 'Precision_AB', 'Precision_BB', 'Recall_AA_tested', 'Recall_AB_tested', 'Recall_BB_tested', 'Number_AA_all', 'Number_AB_all', 'Number_BB_all', 'Recall_AA_all', 'Recall_AB_all', 'Recall_BB_all', 'group']
    performance.to_csv('performance/combined/performance_one_method_%s%s.txt' % (sample, suffix), sep='\t', index = False)


    # Explore the directions
    discrepancy = discrepancy.copy()

    discrepancy_arr = np.array(discrepancy.apply(lambda x: '_'.join((str(x['Originally_called']), str(x['Imputed']))), axis=1))

    AA_BB_to_AB = discrepancy.iloc[np.where(['_1.0' in x for x in discrepancy_arr])[0]]
    AA_BB_to_AB_precision_imputed  = np.sum(AA_BB_to_AB['Imputed'] == AA_BB_to_AB['True']) / float(len(AA_BB_to_AB))
    AA_BB_to_AB_precision_original = np.sum(AA_BB_to_AB['Originally_called'] == AA_BB_to_AB['True']) / float(len(AA_BB_to_AB))

    AB_to_AA_BB = discrepancy.iloc[np.where(['1.0_' in x for x in discrepancy_arr])[0]]
    AB_to_AA_BB_precision_imputed  = np.sum(AB_to_AA_BB['Imputed'] == AB_to_AA_BB['True']) / float(len(AB_to_AA_BB))
    AB_to_AA_BB_precision_original = np.sum(AB_to_AA_BB['Originally_called'] == AB_to_AA_BB['True']) / float(len(AB_to_AA_BB))

    AA_to_BB = discrepancy.iloc[np.where(['1.0' not in x for x in discrepancy_arr])[0]]
    AA_to_BB_precision_imputed  = np.sum(AA_to_BB['Imputed'] == AA_to_BB['True']) / float(len(AA_to_BB))
    AA_to_BB_precision_original = np.sum(AA_to_BB['Originally_called'] == AA_to_BB['True']) / float(len(AA_to_BB))

    stats = pd.DataFrame([sample, len(WGS_df), len(combined), N_notExist_in_WGS, len(only_in_orignal), len(called_in_both) - len(discrepancy), len(discrepancy), len(only_in_imputed)]).transpose()
    stats.columns = ["sample", "Variants_in_WGS","called_and_in_WGS", "called(imputed)_not_in_WGS","only_in_orignal", "called_in_both_consistent", "called_in_both_inconsistently", "only_in_imputed"]

    print('Evalutae genotype imputed - captured in both methods but with disconcordant genotype')
    discrepancy_metrics = pd.DataFrame([sample, len(AA_BB_to_AB), AA_BB_to_AB_precision_imputed, AA_BB_to_AB_precision_original, len(AB_to_AA_BB), AB_to_AA_BB_precision_imputed, AB_to_AA_BB_precision_original, len(AA_to_BB), AA_to_BB_precision_imputed, AA_to_BB_precision_original] + total_number_array).transpose()
    discrepancy_metrics.columns = ["sample", "AA_BB_to_AB", "AA_BB_to_AB_precision_imputed", "AA_BB_to_AB_precision_original", "AB_to_AA_BB", "AB_to_AA_BB_precision_imputed", "AB_to_AA_BB_precision_original", "AA_to_BB", "AA_to_BB_precision_imputed", "AA_to_BB_precision_original", "Number_AA_all", 'Number_AB_all', 'Number_BB_all']

    stats.to_csv('performance/combined/called_in_both_%s%s.txt' % (sample, suffix), sep='\t', index = False)
    discrepancy_metrics.to_csv('performance/combined/called_inconsistently_%s%s.txt' % (sample, suffix), sep='\t', index = False)


if __name__ == '__main__':
    minDP = 2
    sample = sys.argv[1]
    if not os.path.exists('performance/combined'):
        os.makedirs('performance/combined')

    compute_genotype_metrics_called_and_imputed(sample, restrict_to_SNP = False)
    compute_genotype_metrics_called_and_imputed(sample, restrict_to_SNP = True)
