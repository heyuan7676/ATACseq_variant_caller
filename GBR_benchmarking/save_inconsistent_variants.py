
# coding: utf-8

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
from collections import Counter


from utils import *


# In[2]:


def collect_inconsistent_variants(minDP, sample):     
    
    # golden standard
    golden_standard_genotype = readin_golden_standard_genotype(golden_standard_dir, sample)
    
    # genotype caller
    [called_genotype, variants_count] = readin_variant_caller_genotype(VCF_dir, minDP, sample)  
    
    # imputation
    [imputed_genotype, variants_count_imputed] = readin_imputation_genotype(Imputation_dir, minDP, sample)  
    
    # combined the two sources
    called_in_both = called_genotype.merge(imputed_genotype, left_index=True, right_index=True)
    called_both_inconsistent_dat = called_in_both[called_in_both['GT'] != called_in_both['Imputed_GT']]
    called_both_inconsistent_dat = called_both_inconsistent_dat.merge(golden_standard_genotype, left_index = True, right_index=True, how = 'left')
    
    called_both_inconsistent_dat.to_csv('%s/input/%s_minDP%d_inconsistent_variants.txt' % (training_model_dir, sample, minDP), sep='\t')
    
    


# In[ ]:


if __name__ == '__main__':

    minDP = 2
    sample = 'HG00096'

    golden_standard_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/sample_maf_0.05'
    VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/'
    Imputation_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation'

    training_model_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model'
    save_dir = 'Performance'
    
    # golden standard    
    collect_inconsistent_variants(minDP, sample)

