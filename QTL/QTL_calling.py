import sys
import os
import pdb
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from statsmodels.regression import linear_model as sm


def read_in_peaks():
    peak_dat = pd.read_csv('%s/peak_by_sample_matrix_corrected_chr%d.txt' % (PEAK_DIR, CHROMOSOME),
                          sep='\t')
    samples = [x for x in peak_dat.columns if x.startswith('HG')]
    return [peak_dat, samples]

def convert_gt_to_number(arri):
    '''
    convert genotype type A/A to 0, 1, 2 ...
    '''
    arri_idx = np.where(arri!='0')[0]
    nts = np.unique([a for b in [x.split('/') for x in arri[arri_idx]] for a in b])
    # remove the sites with more than two alleles
    if len(nts) > 2:
        return 
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
    '''
    Very slow
    '''

    gt_numerical_dat = []
    print('In total %d rows' % len(gt_dat))
    for i in range(len(gt_dat)):
        arri = np.array(gt_dat[samples].iloc[i])
        token = convert_gt_to_number(arri)
        if token is None:
            continue
        gt_numerical_dat.append(token)
        if i % 500 == 0:
            print('%d rows finished' % i)

    gt_numerical_dat = pd.DataFrame(gt_numerical_dat)
    gt_numerical_dat.columns = samples
    gt_numerical_dat['CHR_POS'] = gt_dat['CHR_POS']
    gt_numerical_dat['CHR'] = gt_dat['CHR']
    gt_numerical_dat['POS'] = gt_dat['POS']
    gt_numerical_dat = gt_numerical_dat[['CHR_POS', 'CHR', 'POS']  + samples]
    
    gt_numerical_dat.to_csv('%s/gt_by_sample_matrix_chr%d_numericalGT.txt' % (Genotype_dir, CHROMOSOME), sep=' ', index = False)
    
    return gt_numerical_dat
    



def readin_genotype():
    ## Read in Genotpye
    gt_dat = pd.read_csv('%s/gt_by_sample_matrix_chr%d.txt' % (Genotype_dir, CHROMOSOME), sep=' ')
    gt_dat = gt_dat[[x for x in gt_dat.columns if 'Unnamed' not in x]]
    samples = [x for x in gt_dat.columns if x.startswith('HG')]

    # remove rows with less than 3 samples
    valid_snps = np.where(np.sum(np.array(gt_dat[samples]) != '0', axis=1) >= 3)[0]
    gt_dat = gt_dat.iloc[valid_snps].reset_index(drop=True)
    
    # remove rows with only one genotype
    validQTLsnps = np.where([len(set(x))>2 for x in np.array(gt_dat[samples])])[0]
    gt_dat = gt_dat.iloc[validQTLsnps].reset_index(drop=True)

    try:
        gt_numerical_dat = pd.read_csv('%s/gt_by_sample_matrix_chr%d_numericalGT.txt' % (Genotype_dir, CHROMOSOME), sep=' ')
    except:
        print('Convert genotype to numerical values')
        gt_numerical_dat = obtain_numerical_gt(gt_dat, samples)
        
    return gt_numerical_dat


def achieve_ll(rowi):
    '''
    Map PL scores back to likelihood
    '''
    
    pl_scores = [list(map(int, x.split(','))) for x in rowi]
    pl_scores = [np.array(x) for x in pl_scores]
    pl_scores = [x[x<10] for x in pl_scores]

    post_pp = [list(map(lambda x: 1/10 ** x, pl_s)) for pl_s in pl_scores]    
    max_post_pp = [np.max(x/np.sum(x)) for x in post_pp]
    
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
    
    post_pp_dat.to_csv('%s/gt_infoLL_by_sample_matrix_chr%d.txt' % (VCF_dir, CHROMOSOME), sep=' ', index = False)
    
    return post_pp_dat




def readin_genotype_info(gt_dat):
    ## Read in Genotpye INFO
    info_dat = pd.read_csv('%s/gt_info_by_sample_matrix_chr%d.txt' % (VCF_dir, CHROMOSOME), sep=' ')
    info_dat = info_dat[[x for x in info_dat.columns if 'Unnamed' not in x]]
    info_dat = info_dat.set_index('CHR_POS').loc[gt_dat['CHR_POS']].reset_index()
    
    samples = [x for x in info_dat.columns if x.startswith('HG')]
    
    try:
        post_pp_dat = pd.read_csv('%s/gt_infoLL_by_sample_matrix_chr%d.txt' % (VCF_dir, CHROMOSOME), sep=' ')
    except:
        print('Derive posterior probability for the genotypes')
        post_pp_dat = derive_ll(info_dat=info_dat, samples=samples)
        
    return post_pp_dat


#### Compute QTL

# Genotype data: gt_numerical_dat
# Peak data: RPKM_dat
# Genotype data weights: post_pp_dat


def compute_QTL_gti_peaki(gt_sample, peak_sample, weight_sample):
    #gt_sample = np.array(gt_numerical_dat[samples].iloc[0])
    #peak_sample = np.array(RPKM_dat[samples].iloc[0])
    #weight_sample = np.array(post_pp_dat[samples].iloc[0])
    valid_samples = np.where(gt_sample!= -1)[0]
    
    y = np.array(peak_sample[valid_samples])
    x = np.array(gt_sample[valid_samples])
    x_weights = np.array(weight_sample[valid_samples])
    
    try:
        wls_model = sm.WLS(y, x, weights = x_weights)
        results = wls_model.fit()
        return results.f_pvalue
    except:
        pdb.set_trace()



def compute_QTL_peaki(peaki, gt_numerical_dat, post_pp_dat):
    QTL_result_peaki = []
    #peaki = RPKM_dat[['START', 'END']].iloc[i]
    [start, end] = [peaki['START'], peaki['END']]
    
    SNPs_close = (np.array(gt_numerical_dat['POS']) < end + WINDOW) & (np.array(gt_numerical_dat['POS']) > start - WINDOW)
    SNPs_close = np.where(SNPs_close)[0]
    
    SNPs_gt = gt_numerical_dat.iloc[SNPs_close]
    SNPs_gt_ll = post_pp_dat.iloc[SNPs_close]
    
    for gti in range(len(SNPs_gt)):
        print(peaki['PEAK'], SNPs_gt.iloc[gti]['CHR_POS'])
        pval = compute_QTL_gti_peaki(SNPs_gt[samples].iloc[gti], peaki[samples], SNPs_gt_ll[samples].iloc[gti])
        if(np.isnan(pval)):
            break
        snpid = SNPs_gt.iloc[gti]['CHR_POS']
        peakid = peaki['PEAK']
        QTL_result_peaki.append([peakid, snpid, pval])
        
    return QTL_result_peaki



def compute_QTLs(WINDOW, peak_df, genotype_df, weight_df):

    QTL_results = []
    for p in range(len(peak_df)):
        token = compute_QTL_peaki(peak_df.iloc[p], genotype_df, weight_df)
        QTL_results.append([a for b in token for a in b])
        if p % 1000 == 0:
            print('%d peaks finished' % p)
        
    QTL_results = [a for b in QTL_results for a in b]
    QTL_results = np.reshape(np.array(QTL_results), [int(len(QTL_results)/3), 3])
    QTL_dat = pd.DataFrame(QTL_results)
    QTL_dat.columns = ['PeakID', 'CHR_POS', 'p-value']
    QTL_dat['p-value'] = list(map(float, QTL_dat['p-value']))
   
    print('Tested %d associations for chromosome%d' % (len(QTL_dat), CHROMOSOME))
    print("") 
    QTL_dat.to_csv('%s/CHR%d_acQTLs_WINDOW_%dkb.txt' % (QTL_dir, CHROMOSOME, WINDOW/1000), sep='\t', index = False)


if __name__ == "__main__":
    CHROMOSOME = int(sys.argv[1])
    WINDOW = int(sys.argv[2])

    PEAK_DIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks'
    BAM_DIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/first_pass_bqsr'
    Genotype_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Called_GT/'
    VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files'
    QTL_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/QTLs'

    [PEAK_DAT, samples] = read_in_peaks()
    GT_DAT = readin_genotype()
    WEIGHT_DAT = readin_genotype_info(gt_dat=GT_DAT)

    compute_QTLs(WINDOW, PEAK_DAT, GT_DAT, WEIGHT_DAT)

