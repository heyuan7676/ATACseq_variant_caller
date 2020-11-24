import sys
import time
import os
import pdb
import numpy as np
import pandas as pd


def convert_gt_to_number(arri):
    '''
    convert genotype type A/A to 0, 1, 2 ...
    '''
    arri_idx = np.where(arri!='0')[0]
    nts = np.unique([a for b in [x.split('/') for x in arri[arri_idx]] for a in b])
    # remove the sites with more than two alleles
    if len(nts) > 2:
        return [-1] * len(arri)
    nts_dict = {}
    i = 0
    for ni in nts:
        nts_dict[ni] = i
        i += 1
    numeric_gt = [np.sum([nts_dict[a] for a in x.split('/')]) for x in arri[arri_idx]]
    numeric_gt_arr = np.ones(len(arri)) * (-1)
    numeric_gt_arr[arri_idx] = numeric_gt
    return numeric_gt_arr


def obtain_numerical_gt(gt_dat, samples):

    gt_numerical_dat = []
    print('In total %d variants are called' % len(gt_dat))
    start = time.time()
    gt_numerical_dat = list(map(convert_gt_to_number, np.array(gt_dat[samples])))
    end = time.time()
    print('    Used %f s' % (end - start))
    
    gt_numerical_dat = pd.DataFrame(gt_numerical_dat)

    # remove rows with more than two alleles
    validQTLsnps = np.where([len(set(x))>1 for x in np.array(gt_numerical_dat)])[0]
    gt_numerical_dat = gt_numerical_dat.iloc[validQTLsnps].reset_index(drop=True)
    gt_dat = gt_dat.iloc[validQTLsnps].reset_index(drop=True)
 
    gt_numerical_dat.columns = samples
    gt_numerical_dat['CHR_POS'] = gt_dat['CHR_POS']
    gt_numerical_dat['CHR'] = gt_dat['CHR']
    gt_numerical_dat['POS'] = gt_dat['POS']
    gt_numerical_dat = gt_numerical_dat[['CHR_POS', 'CHR', 'POS']  + list(samples)]
    
    return gt_numerical_dat
    


def readin_genotype(Genotype_dir, chromosome, samples):
    ## Read in Genotpye
    gt_dat = pd.read_csv('%s/gt_by_sample_matrix_chr%d.txt' % (Genotype_dir, chromosome), sep=' ')
    gt_dat = gt_dat[[x for x in gt_dat.columns if 'Unnamed' not in x]]
    gt_dat = gt_dat.replace('./.', '0')
    gt_dat = gt_dat.replace(0, '0')
    
    # match sample order in genotype matrix with the peak matrix
    gt_dat = gt_dat[['CHR_POS', 'CHR', 'POS'] + samples]

    # remove rows with less than 3 samples
    valid_snps = np.where(np.sum(np.array(gt_dat[samples]) != '0', axis=1) >= 3)[0]
    gt_dat = gt_dat.iloc[valid_snps].reset_index(drop=True)
    
    # remove rows with only one genotype
    validQTLsnps = np.where([len(set(x))>2 for x in np.array(gt_dat[samples])])[0]
    gt_dat = gt_dat.iloc[validQTLsnps].reset_index(drop=True)

    print('Convert genotype to numerical values')
    gt_numerical_dat = obtain_numerical_gt(gt_dat, samples)

    # again remove rows with only one genotype (ie. A/T and T/A)
    validQTLsnps = np.where([len(set(x))>2 for x in np.array(gt_numerical_dat[samples])])[0]
    gt_numerical_dat = gt_numerical_dat.iloc[validQTLsnps].reset_index(drop=True)
        
    return gt_numerical_dat


def achieve_ll(rowi):
    '''
    Map PL scores back to likelihood
    '''
    pl_scores = [list(map(int, str(x).split(','))) for x in rowi]
    pl_scores = [np.array(x) for x in pl_scores]
    pl_scores = [x[x<10] for x in pl_scores]

    post_pp = [list(map(lambda x: 1/10 ** x, pl_s)) for pl_s in pl_scores]    
    max_post_pp = [np.max(np.array(x)/float(np.sum(x))) for x in post_pp]
    
    return max_post_pp


def derive_ll(info_dat, samples):
    info_arr = np.array(info_dat[samples])
    idx_score = np.where(info_arr!='0')
    post_pp_arr = achieve_ll(info_arr[idx_score])

    post_pp = np.zeros(info_arr.shape)
    post_pp[idx_score] = post_pp_arr
    
    post_pp_dat = pd.DataFrame(post_pp)
    post_pp_dat.columns = samples
    post_pp_dat['CHR_POS'] = info_dat['CHR_POS']
    post_pp_dat['CHR'] = info_dat['CHR']
    post_pp_dat['POS'] = info_dat['POS']
    post_pp_dat = post_pp_dat[['CHR_POS', 'CHR', 'POS']  + samples]
    
    return post_pp_dat




def readin_genotype_info(gt_dat, VCF_dir, chromosome, samples):
    ## Read in Genotpye INFO
    info_dat = pd.read_csv('%s/gt_info_by_sample_matrix_chr%d.txt' % (VCF_dir, chromosome), sep=' ')
    info_dat = info_dat[[x for x in info_dat.columns if 'Unnamed' not in x]]
    info_dat = info_dat.set_index('CHR_POS').loc[gt_dat['CHR_POS']].reset_index()
    info_dat = info_dat.replace(0, '0')
   
    # match sample order in genotype matrix with the peak matrix
    info_dat = info_dat[['CHR_POS', 'CHR', 'POS'] + samples]
 
    print('Derive posterior probability for the genotypes')
    start = time.time()
    post_pp_dat = derive_ll(info_dat=info_dat, samples=samples)
    end = time.time()
    print('    Used %f seconds' % (end - start))
        
    return post_pp_dat


