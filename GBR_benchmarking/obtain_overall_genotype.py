import os
import sys
import numpy as np
import pandas as pd

from utils import *
from scipy.stats import spearmanr



def compute_metrics(df_two_col):
    df_two_col = df_two_col[~df_two_col['True_GT'].isnull()]
    spearman_cor_i = spearmanr(df_two_col['Dosage'], df_two_col['True_GT'])[0]
    mse_i = np.mean([x**2 for x in df_two_col['Dosage'] - df_two_col['True_GT']])
    rowi = [spearman_cor_i, mse_i]  
    for gt in [0, 1, 2]:
        df_two_col_gt = df_two_col[df_two_col['True_GT'] == gt]
        mse_i = np.mean([x**2 for x in df_two_col_gt['Dosage'] - df_two_col_gt['True_GT']])
        rowi.append(mse_i)
    rowi.append(len(df_two_col))
    return rowi




def evaluate_variants_all(minDP, sample):

    performance = []

    # predicted
    predicted_dat = pd.read_csv('%s/output/%s_minDP%d_variants_THR%.1f.txt' % (training_model_dir, sample, minDP, THR), sep='\t', index_col = 0)

    # genotype caller
    [called_genotype, variants_count] = readin_variant_caller_genotype(VCF_dir, minDP, sample)
    called_genotype['True_GT'] = golden_standard_genotype.loc[np.array(called_genotype.index)]['True_GT']

    rowi = [[sample, minDP, 'GC']] 
    rowi.append(compute_metrics(called_genotype))
    performance.append([a for b in rowi for a in b])

    # imputation
    [imputed_genotype, variants_count_imputed] = readin_imputation_genotype(Impute_dir, minDP, sample)
    imputed_genotype['True_GT'] = golden_standard_genotype.loc[np.array(imputed_genotype.index)]['True_GT']
    imputed_genotype = imputed_genotype[['Imputed_dosage', 'True_GT']]
    imputed_genotype.columns = ['Dosage', 'True_GT']

    rowi = [[sample, minDP, 'Imputation']] 
    rowi.append(compute_metrics(imputed_genotype))
    performance.append([a for b in rowi for a in b])

    # combined
    set1 = np.array(list(set(called_genotype.index) - set(predicted_dat.index)))
    set2 = np.array(list(set(imputed_genotype.index) - set(predicted_dat.index)))

    called_genotype = called_genotype.loc[set1]
    imputed_genotype = imputed_genotype.loc[set2]

    for col in ['Imputed_dosage', 'Dosage', 'Y_random_forest', 'Y_linear_regression', 'Y_logistic_regression', 'Y_ordinal_regression']:

        rowi = [[sample, minDP, col]]
        pred_dat = predicted_dat[[col, 'True_GT']]
        pred_dat.columns = ['Dosage', 'True_GT']
        combined = called_genotype.append(imputed_genotype).append(pred_dat)
        #combined.to_csv('%s/%s_minDP%d_variants_THR%.1f_%s.txt' % (integration_dir, sample, minDP, THR, col), sep='\t')

        rowi.append(compute_metrics(combined))
        performance.append([a for b in rowi for a in b])

    performance = pd.DataFrame(performance)
    performance.columns = ['sample', 'minDP', 'method', 'correlation', 'MSE', 'MSE_AA', 'MSE_AB', 'MSE_BB', 'Number']

    return performance




 
if __name__ == '__main__':

    sample = sys.argv[1]
    THR = 0.0

    save_dir = 'Performance'
    save_dat = pd.DataFrame()

    golden_standard_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'
    golden_standard_genotype = readin_golden_standard_genotype(golden_standard_dir, sample)
    if golden_standard_genotype is None:
        sys.exit()


    for minDP in [2,3,4,5,6,8,10]:
        training_model_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model/minDP%d' % minDP
        integration_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Integration/minDP%d' % minDP
        VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/minDP%d' % minDP
        Impute_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation/minDP%d' % minDP
        suffix = '.filtered.minDP%d.imputed.dosage_genotype.bed' % minDP

        os.makedirs(integration_dir, exist_ok = True)
        try:
            evaluate_variants_all_minDP = evaluate_variants_all(minDP, sample)
            save_dat = save_dat.append(evaluate_variants_all_minDP)
        except:
            print('%s minDP%d prediction file does not exist' % (sample, minDP))
            continue

    save_dat.to_csv('%s/Overall_performance_%s.txt' % (save_dir, sample), sep='\t', index = False)




