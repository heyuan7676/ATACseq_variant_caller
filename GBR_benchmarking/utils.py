import os
import sys
import numpy as np
import pandas as pd
from collections import Counter
import pdb


def readin_golden_standard_genotype(golden_standard_dir, sample):    
    print('Read in golden standard variants genotype data')
    chromosome = 22
    golden_standard_fn = '%s/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.vcf' % (golden_standard_dir, chromosome)
    for l in open(golden_standard_fn, 'r').readlines():
        if 'CHROM' in l:
            break
    try:
        sample_col = l.rstrip().split('\t').index(sample)
    except:
        print('%s not present in phase3 data from 1000 Genome' % sample)
        return

    golden_standard_variants = []
    golden_standard_gt = []
    for chromosome in range(1, 23):
        golden_standard_fn = '%s/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.vcf' % (golden_standard_dir, chromosome)
        golden_standard_genotype = pd.read_csv(golden_standard_fn, comment="#", sep='\t', header = None, usecols=[0, 1, sample_col])
        golden_standard_genotype.columns = ['CHR', 'POS', 'GT']
	multiallelic_variants = np.array( golden_standard_genotype['POS'][golden_standard_genotype['POS'].duplicated()])
	biallelic_variants = np.array(list(set(golden_standard_genotype['POS']) - set(multiallelic_variants)))
	golden_standard_genotype = golden_standard_genotype.set_index('POS').loc[biallelic_variants].reset_index()
        golden_standard_gt.append([np.sum(list(map(int, x.split('|')))) for x in golden_standard_genotype['GT']])
        golden_standard_variants.append(['_'.join(map(str, x)) for x in zip(np.array(golden_standard_genotype['CHR']), np.array(golden_standard_genotype['POS']))])
    
    golden_standard_gt = [a for b in golden_standard_gt for a in b]
    golden_standard_variants = [a for b in golden_standard_variants for a in b]
    
    golden_standard_dat = pd.DataFrame({"True_GT": golden_standard_gt})
    golden_standard_dat.index = golden_standard_variants
        
    return golden_standard_dat


# In[3]:


def readin_variant_caller_genotype(VCF_dir, minDP, sample):    
    print('Read in variants with minDP >= %d for %s ...' % (minDP, sample))
    genotype_dosage_filename = '%s/minDP%d/%s.filtered.minDP%d.recode.dosage_genotype.bed' % (VCF_dir, minDP, sample, minDP)
    genotype_dosage = pd.read_csv(genotype_dosage_filename, sep='\t', index_col = 0)    
    genotype_dosage = genotype_dosage[genotype_dosage['Dosage'] != -1]
    print('Done, read %d bi-allelic variants' % len(genotype_dosage))
    
    called_variants = pd.DataFrame.from_dict(Counter(genotype_dosage['GT']), orient='index')
    
    return [genotype_dosage, called_variants]


# In[195]:


def readin_imputation_genotype(Imputation_dir, minDP,sample): 
    print('Read in imputed variants with minDP >= %d for %s' % (minDP, sample))
    imputation_filename = '%s/minDP%d/%s.filtered.minDP%d.imputed.dosage_genotype.bed' % (Imputation_dir, minDP, sample, minDP)
    genotype_imputed = pd.read_csv(imputation_filename, sep='\t', index_col = 0, header = None)
    genotype_imputed.columns = ['GT', 'Imputed_dosage']
    genotype_imputed['Imputed_GT'] = [np.sum(list(map(int, x.split('|')))) for x in np.array(genotype_imputed['GT'])]
    genotype_imputed = genotype_imputed[['Imputed_dosage','Imputed_GT']]
    genotype_imputed = genotype_imputed[~genotype_imputed.index.duplicated()]
    print('Done, read %d variants' % len(genotype_imputed))
    
    called_variants = pd.DataFrame.from_dict(Counter(genotype_imputed['Imputed_GT']), orient='index')
    
    return [genotype_imputed, called_variants]



def obtain_confustion_matrix_r2(combined_dat, golden_variants, called_variants, group):
    combined_dat = combined_dat.copy()
    combined_dat.columns = ['Dosage', 'GT', 'True_GT']
        
    correct_dat = combined_dat[combined_dat['GT'] == combined_dat['True_GT']]
    correct_variants = pd.DataFrame.from_dict(Counter(correct_dat['GT']), orient='index')
    
    confusion_matrix = golden_variants.merge(called_variants, left_index=True, right_index=True).merge(correct_variants, left_index=True, right_index=True)
    confusion_matrix.columns = ['Golden_standard', 'Called', 'Called_and_True']    
    
    confusion_matrix['Precision'] = confusion_matrix['Called_and_True'] / confusion_matrix['Called']
    confusion_matrix['Recall'] = confusion_matrix['Called_and_True'] / confusion_matrix['Golden_standard']
    
    confusion_matrix['True_GT'] = list(map(int, confusion_matrix.index))
    
    combined_dat = combined_dat[~combined_dat['True_GT'].isnull()]
    pearson_correlation = np.corrcoef(np.array(combined_dat['Dosage']), 
                                      np.array(combined_dat['True_GT']))[0][1]
    confusion_matrix['Pearson_correlation'] = pearson_correlation
    confusion_matrix['Group'] = group

    return confusion_matrix
