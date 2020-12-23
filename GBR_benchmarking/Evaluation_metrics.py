
# coding: utf-8

# In[5]:


import os
import sys
import numpy as np
import pandas as pd
from collections import Counter

from utils import *

from matplotlib import pyplot as plt


# In[2]:


def compare_called_variants_to_golden_standard(golden_standard_genotype, minDP, sample, plot_dosage = False):      
    
    # golden standard
    golden_variants = pd.DataFrame.from_dict(Counter(golden_standard_genotype['True_GT']), orient='index')
    
    # genotype caller
    [called_genotype, variants_count] = readin_variant_caller_genotype(VCF_dir, minDP, sample)  
    called_genotype_combined = called_genotype.merge(golden_standard_genotype, left_index = True, right_index=True)
    conf_called = obtain_confustion_matrix_r2(called_genotype_combined, golden_variants, variants_count)
    conf_called.to_csv('%s/%s_minDP%d_genotype_caller_number_precision_recall_r.txt' % (save_dir, sample, minDP), index = False)
    
    # imputation
    [imputed_genotype, variants_count_imputed] = readin_imputation_genotype(Imputation_dir, minDP, sample)   
    imputed_genotype_combined = imputed_genotype.merge(golden_standard_genotype, left_index = True, right_index=True)
    conf_imputed = obtain_confustion_matrix_r2(imputed_genotype_combined, golden_variants, variants_count_imputed)
    conf_imputed.to_csv('%s/%s_minDP%d_imputation_number_precision_recall_r.txt' % (save_dir, sample, minDP), index = False)
    
    if plot_dosage:
        plt.figure()
        plt.hist(called_genotype['Dosage'], histtype = 'step', label = 'Genotype caller', bins = 50, linewidth = 2, color = 'orange')
        plt.hist(imputed_genotype['Imputed_dosage'], histtype = 'step', label = 'Imputation', bins = 50, linewidth = 2, color = 'blue')
        plt.xlabel('Dosage')
        plt.ylabel('Number of variants called')
        plt.title('Distribution of dosage for %s' % sample)
        plt.legend()
        plt.xlim([-0.1, 2.1])
        plt.tight_layout()
        plt.savefig('%s/%s_minDP%d_dosage_distribution.png' % (save_dir, sample, minDP))
        plt.close()    
        
    # combined the two sources
    called_in_both = called_genotype.merge(imputed_genotype, left_index=True, right_index=True)
    called_only_from_gc = list(set(called_genotype.index) - set(imputed_genotype.index))
    called_only_from_imputation = list(set(imputed_genotype.index) - set(called_genotype.index))
    called_both_consistent = called_in_both[called_in_both['GT'] == called_in_both['Imputed_GT']]
    called_both_inconsistent = called_in_both[called_in_both['GT'] != called_in_both['Imputed_GT']]
    
    variants_compostion = [len(called_only_from_gc), len(called_only_from_imputation), 
                           len(called_both_consistent), len(called_both_inconsistent)]
    variants_compostion = pd.DataFrame(variants_compostion).transpose()
    variants_compostion.columns = ['Only_from_gc', 'Only_from_imputation', 'Two_consistent', 'Two_inconsistent']
    variants_compostion.to_csv('%s/%s_minDP%d_variants_numbers.txt' % (save_dir, sample, minDP), sep='\t', index = False)

    dat1 = called_genotype.loc[np.array(called_only_from_gc)]
    dat2 = imputed_genotype.loc[np.array(called_only_from_imputation)]
    
    dat31 = called_genotype.loc[np.array(called_both_consistent.index)]
    dat32 = imputed_genotype.loc[np.array(called_both_consistent.index)]
    dat3 = pd.DataFrame({"Dosage": (dat31['Dosage'] + dat32['Imputed_dosage']) / 2, "GT": dat32['Imputed_GT']}) 

    dat1.columns = ['Dosage', 'GT']
    dat2.columns = ['Dosage', 'GT']
    dat3.columns = ['Dosage', 'GT']

    # 1). discard all inconsistent calls
    intergrated_dat = dat1.append(dat2).append(dat3)
    intergrated_variants = pd.DataFrame.from_dict(Counter(intergrated_dat['GT']), orient='index')
    
    intergrated_genotype_combined = intergrated_dat.merge(golden_standard_genotype, left_index = True, right_index=True)
    conf_overall_no_inconsistent = obtain_confustion_matrix_r2(intergrated_genotype_combined, golden_variants, intergrated_variants)
    conf_overall_no_inconsistent.to_csv('%s/%s_minDP%d_combined_no_inconsistent_number_precision_recall_r.txt' % (save_dir, sample, minDP), index = False)
    
    # 2). add variant learned from the model
    try:
        dat4 = pd.read_csv('%s/output/%s_minDP%d_inconsistent_variants.txt' % (training_model_dir, sample, minDP), sep='\t', index_col = 0)
    except:
        return
    
    dat4.columns = ['Dosage', 'GT']
    
    intergrated_dat = intergrated_dat.append(dat4)
    intergrated_variants = pd.DataFrame.from_dict(Counter(intergrated_dat['GT']), orient='index')
    
    intergrated_genotype_combined = intergrated_dat.merge(golden_standard_genotype, left_index = True, right_index=True)
    conf_overall_with_inconsistent = obtain_confustion_matrix_r2(intergrated_genotype_combined, golden_variants, intergrated_variants)
    conf_overall_with_inconsistent.to_csv('%s/%s_minDP%d_combined_with_inconsistent_number_precision_recall_r.txt' % (save_dir, sample, minDP), index = False)
        
    return
    


# In[6]:


if __name__ == '__main__':

    minDP = 2
    sample = 'HG00096'
    plot_dosage = True

    golden_standard_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/sample_maf_0.05'
    VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/'
    Imputation_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation'

    training_model_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model'
    save_dir = 'Performance'
    
    # golden standard    
    golden_standard_genotype = readin_golden_standard_genotype(golden_standard_dir, sample)
    compare_called_variants_to_golden_standard(golden_standard_genotype, minDP, sample, plot_dosage = plot_dosage)

