
# coding: utf-8

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
from collections import Counter


# In[ ]:


golden_standard_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'

def readin_golden_standard_genotype(sample):    
    print('Read in golden standard variants genotype data')
    chromosome = 22
    golden_standard_fn = '%s/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.vcf' % (golden_standard_dir, chromosome)
    for l in open(golden_standard_fn, 'r').readlines():
        if 'CHROM' in l:
            break
    sample_col = l.rstrip().split('\t').index(sample)
    
    golden_standard_variants = []
    golden_standard_gt = []
    for chromosome in range(1, 23):
        golden_standard_fn = '%s/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.vcf' % (golden_standard_dir, chromosome)
        golden_standard_genotype = pd.read_csv(golden_standard_fn, comment="#", sep='\t', header = None, usecols=[0, 1, sample_col])
        golden_standard_genotype.columns = ['CHR', 'POS', 'GT']
        golden_standard_gt.append([np.sum(list(map(int, x.split('|')))) for x in golden_standard_genotype['GT']])
        golden_standard_variants.append(['_'.join(map(str, x)) for x in zip(np.array(golden_standard_genotype['CHR']), np.array(golden_standard_genotype['POS']))])
    
    golden_standard_gt = [a for b in golden_standard_gt for a in b]
    golden_standard_variants = [a for b in golden_standard_variants for a in b]
    
    golden_standard_dat = pd.DataFrame({"True_GT": golden_standard_gt})
    golden_standard_dat.index = golden_standard_variants
        
    return golden_standard_dat


# In[6]:


VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/'

def readin_variant_caller_genotype(minDP, sample):    
    print('Read in variants with minDP >= %d for %s' % (minDP, sample))
    genotype_dosage_filename = '%s/minDP%d/%s.filtered.minDP%d.recode.dosage_genotype.bed' % (VCF_dir, minDP, sample, minDP)
    genotype_dosage = pd.read_csv(genotype_dosage_filename, sep='\t', index_col = 0)    
    return genotype_dosage


# In[7]:


Imputation_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation'

def readin_imputation_genotype(minDP,sample): 
    imputation_filename = '%s/minDP%d/%s.filtered.minDP%d.imputed.dosage_genotype.bed' % (Imputation_dir, minDP, sample, minDP)
    genotype_imputed = pd.read_csv(imputation_filename, sep=' ', index_col = 0, header = None)
    genotype_imputed.columns = ['GT', 'Imputed_dosage']
    genotype_imputed['Imputed_GT'] = [np.sum(list(map(int, x.split('|')))) for x in np.array(genotype_imputed['GT'])]
    genotype_imputed = genotype_imputed[['Imputed_GT', 'Imputed_dosage']]
    
    return genotype_imputed


# In[11]:


def obtain_confustion_matrix_r2(combined_dat):
    recovered_dat = combined_dat[~combined_dat['GT'].isnull()]
    recovered_dat = recovered_dat[combined_dat['Dosage'] != -1]
    called_variants = pd.DataFrame.from_dict(Counter(recovered_dat['GT']), orient='index')
    pearson_correlation = np.corrcoef(np.array(recovered_dat['Dosage']), 
                                      np.array(recovered_dat['True_GT']))[1][1]
    
    correct_dat = recovered_dat[recovered_dat['GT'] == recovered_dat['True_GT']]
    correct_variants = pd.DataFrame.from_dict(Counter(correct_dat['GT']), orient='index')
    
    confusion_matrix = golden_variants.merge(called_variants, left_index=True, right_index=True).merge(correct_variants, left_index=True, right_index=True)
    confusion_matrix.columns = ['Golden_standard', 'Called', 'Called_and_True']    
    
    confusion_matrix['Precision'] = confusion_matrix['Called_and_True'] / confusion_matrix['Called']
    confusion_matrix['Recall'] = confusion_matrix['Called_and_True'] / confusion_matrix['Golden_standard']

    return [confusion_matrix, pearson_correlation]    


# In[ ]:


minDP = 2
sample = 'HG00096'


# In[7]:


#def compare_called_variants_to_golden_standard(minDP, sample):  

if 1:
    
    # golden standard
    golden_standard_genotype = readin_golden_standard_genotype(sample = sample)
    golden_variants = pd.DataFrame.from_dict(Counter(golden_standard_genotype['True_GT']), orient='index')
    
    # genotype caller
    called_genotype = readin_variant_caller_genotype(minDP=minDP, sample=sample)  
    called_genotype_combined = called_genotype.merge(golden_standard_genotype, left_index = True, right_index=True, how = 'right')
    [conf_called, correlation_called] = obtain_confustion_matrix_r2(called_genotype_combined)
    
    # imputation
    imputed_genotype = readin_imputation_genotype(minDP=minDP, sample=sample)  
    imputed_genotype_combined = imputed_genotype.merge(golden_standard_genotype, left_index = True, right_index=True, how = 'right')
    [conf_imputed, correlation_imputed] = obtain_confustion_matrix_r2(imputed_genotype_combined)

    


# In[8]:


minDP = 2
sample = 'HG00096'


# In[9]:


called_genotype = readin_variant_caller_genotype(minDP=minDP, sample=sample)  


# In[10]:


called_genotype

