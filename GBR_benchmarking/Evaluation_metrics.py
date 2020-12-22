import os
import sys
import numpy as np
import pandas as pd
from collections import Counter



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



VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/'

def readin_variant_caller_genotype(minDP, sample):    
    print('Read in variants with minDP >= %d for %s ...' % (minDP, sample))
    genotype_dosage_filename = '%s/minDP%d/%s.filtered.minDP%d.recode.dosage_genotype.bed' % (VCF_dir, minDP, sample, minDP)
    genotype_dosage = pd.read_csv(genotype_dosage_filename, sep='\t', index_col = 0)    
    print('Done, read %d variants' % len(genotype_dosage))
    return genotype_dosage




Imputation_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation'

def readin_imputation_genotype(minDP,sample): 
    print('Read in imputed variants with minDP >= %d for %s' % (minDP, sample))
    imputation_filename = '%s/minDP%d/%s.filtered.minDP%d.imputed.dosage_genotype.bed' % (Imputation_dir, minDP, sample, minDP)
    genotype_imputed = pd.read_csv(imputation_filename, sep=' ', index_col = 0, header = None)
    genotype_imputed.columns = ['GT', 'Imputed_dosage']
    genotype_imputed['Imputed_GT'] = [np.sum(list(map(int, x.split('|')))) for x in np.array(genotype_imputed['GT'])]
    genotype_imputed = genotype_imputed[['Imputed_dosage','Imputed_GT']]
    print('Done, read %d variants' % len(genotype_imputed))
    
    return genotype_imputed


# In[6]:


def obtain_confustion_matrix_r2(combined_dat, golden_variants):
    combined_dat = combined_dat.copy()
    combined_dat.columns = ['Dosage', 'GT', 'True_GT']
        
    recovered_dat = combined_dat[~combined_dat['GT'].isnull()]
    recovered_dat = recovered_dat[recovered_dat['Dosage'] != -1]
    called_variants = pd.DataFrame.from_dict(Counter(recovered_dat['GT']), orient='index')
    pearson_correlation = np.corrcoef(np.array(recovered_dat['Dosage']), 
                                      np.array(recovered_dat['True_GT']))[0][1]
    
    correct_dat = recovered_dat[recovered_dat['GT'] == recovered_dat['True_GT']]
    correct_variants = pd.DataFrame.from_dict(Counter(correct_dat['GT']), orient='index')
    
    confusion_matrix = golden_variants.merge(called_variants, left_index=True, right_index=True).merge(correct_variants, left_index=True, right_index=True)
    confusion_matrix.columns = ['Golden_standard', 'Called', 'Called_and_True']    
    
    confusion_matrix['Precision'] = confusion_matrix['Called_and_True'] / confusion_matrix['Called']
    confusion_matrix['Recall'] = confusion_matrix['Called_and_True'] / confusion_matrix['Golden_standard']
    confusion_matrix['Pearson_correlation'] = pearson_correlation
    confusion_matrix['True_GT'] = list(map(int, confusion_matrix.index))

    return confusion_matrix


# In[7]:


save_dir = 'Performance'
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(font_scale = 1.2)
sns.set_style('whitegrid')

def compare_called_variants_to_golden_standard(minDP, sample, plot_dosage = False):      
    # golden standard
    golden_standard_genotype = readin_golden_standard_genotype(sample = sample)
    golden_variants = pd.DataFrame.from_dict(Counter(golden_standard_genotype['True_GT']), orient='index')
    
    # genotype caller
    called_genotype = readin_variant_caller_genotype(minDP=minDP, sample=sample)  
    called_genotype_combined = called_genotype.merge(golden_standard_genotype, left_index = True, right_index=True, how = 'right')
    conf_called = obtain_confustion_matrix_r2(called_genotype_combined, golden_variants)
    conf_called.to_csv('%s/%s_minDP%d_genotype_caller_number_precision_recall_r.txt' % (save_dir, sample, minDP), index = False)
    
    # imputation
    imputed_genotype = readin_imputation_genotype(minDP=minDP, sample=sample)  
    imputed_genotype_combined = imputed_genotype.merge(golden_standard_genotype, left_index = True, right_index=True, how = 'right')
    conf_imputed = obtain_confustion_matrix_r2(imputed_genotype_combined, golden_variants)
    conf_imputed.to_csv('%s/%s_minDP%d_imputation_number_precision_recall_r.txt' % (save_dir, sample, minDP), index = False)
    
    # combined the two sources
    called_genotype_combined = called_genotype_combined[~called_genotype_combined['Dosage'].isnull()]
    called_genotype_combined = called_genotype_combined[called_genotype_combined['Dosage'] != -1]

    imputed_genotype_combined = imputed_genotype_combined[~imputed_genotype_combined['Imputed_dosage'].isnull()]
    imputed_genotype_combined = imputed_genotype_combined[imputed_genotype_combined['Imputed_dosage'] != -1]

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

    called_in_both = called_genotype_combined.merge(imputed_genotype_combined, left_index=True, right_index=True)
    
    return called_in_both


# In[8]:


minDP = 2
sample = 'HG00096'

dat = compare_called_variants_to_golden_standard(minDP = minDP, sample=sample, plot_dosage=True)

