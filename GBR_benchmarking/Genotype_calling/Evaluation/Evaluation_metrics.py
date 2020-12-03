import pdb
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ranksums

def read_in_WGS_GT(sample, assembly = 'GRCh38'):
    print('Read in WGS data...')
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/maf005'

    chromosome = 22
    if assembly == 'GRCh38':
        suffix = 'genotypes.tsv'
    elif assembly == 'GRCh37':
        suffix = 'GRCh37.genotypes.tsv'

    WGS_fn = '%s/1k_genome_chr%d.%s' % (WGS_dir, chromosome, suffix)
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t', nrows = 10)
    try:
        col_to_read = list(WGS_result.columns).index(sample)
    except:
        print('%s does not have genotype data' % sample)
        return 
    
    WGS_result = pd.DataFrame()

    for chromosome in range(21, 23):
        WGS_fn = '%s/1k_genome_chr%d.%s' % (WGS_dir, chromosome, suffix)
        token = pd.read_csv(WGS_fn, 
                            comment = '$', 
                            sep='\t', 
                            usecols=[0,1,2,col_to_read])
        WGS_result = WGS_result.append(token)    
    WGS_result = WGS_result.drop_duplicates() 
    print("done\n")
    return WGS_result


def obtain_atac_variants_df(sample, WGS_result, restrict_to_SNP = True, return_df = True, Imputed = False, minDP = 2, return_metric = True):

    print('Read in genotype data called from ATAC-seq reads...')
    WGS_result = WGS_result.copy()

    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie'
    if Imputed:
        SNP_calling_dir = '%s/Imputation/%s' % (root_dir, sample)
        SNP_called_fn = '%s/%s.imputed.GRCh38.biallelic.genotype.txt' % (SNP_calling_dir, sample)
    else:
        SNP_calling_dir = '%s/Called_GT/minDP%d' % (root_dir, minDP)
        SNP_called_fn = '%s/%s.filtered.genotype.minDP%d.txt' % (SNP_calling_dir, sample, minDP)

    SNP_called = pd.read_csv(SNP_called_fn, comment="$", sep='\t', nrows = 10)
    try:
        col_to_read = np.where([sample in x for x in SNP_called.columns])[0][0]
    except:
        print('%s does not called genotpye' % sample)
        return

    SNP_called = pd.read_csv(SNP_called_fn, comment="$", sep='\t', usecols=[0, 1, 2, col_to_read], low_memory=False)
    SNP_called.columns = list(SNP_called.columns[:3]) + ['%s_called' % sample]
    SNP_called['#CHROM'] = SNP_called['#CHROM'].apply(str)
    WGS_result['#CHROM'] = WGS_result['#CHROM'].apply(str)

    # restrict to SNP
    if restrict_to_SNP:
        WGS_result = WGS_result.iloc[np.where([len(x) == 1 for x in WGS_result['REF']])[0]]
        WGS_result = WGS_result.iloc[np.where([len(x.split('/')[0]) == 1 for x in WGS_result[sample]])[0]]
        WGS_result = WGS_result.iloc[np.where([len(x.split('/')[1]) == 1 for x in WGS_result[sample]])[0]]

        SNP_called = SNP_called.iloc[np.where([len(x.split('/')[1]) == 1 for x in SNP_called['%s_called' % sample]])[0]]
        SNP_called = SNP_called.iloc[np.where([len(x.split('/')[1]) == 1 for x in SNP_called['%s_called' % sample]])[0]]

    # variants with some read
    SNP_called = SNP_called.drop_duplicates()
    intersection_SNPs = WGS_result.merge(SNP_called, on=['#CHROM', 'POS'], how = 'outer')
    #assert len(intersection_SNPs) == len(WGS_result)
    intersection_SNPs.loc[intersection_SNPs['REF_y'].isnull(),'%s_called' % sample] = 'N/N'
    intersection_SNPs.loc[intersection_SNPs['REF_x'].isnull(),'%s' % sample] = 'N/N'

    N1 = sum(~intersection_SNPs['REF_y'].isnull())
    N2 = sum(intersection_SNPs['REF_x'].isnull())
    print("Among %d variants identified by ATAC-seq reads, %d (%.3f) are not included in WGS data" % (N1, N2, float(N2) / N1))

    # match for the HT order
    intersection_SNPs['%s_makeup' % sample] = ['/'.join(x.split('/')[::-1]) for x in intersection_SNPs[sample]]

    HT = np.where(intersection_SNPs[sample] != intersection_SNPs['%s_makeup' % sample])[0]
    match_idx = np.where(intersection_SNPs['%s_called' % sample] == intersection_SNPs['%s_makeup' % sample])[0]
    need_to_flip = np.intersect1d(HT, match_idx)
    intersection_SNPs.loc[need_to_flip,sample] = intersection_SNPs.loc[need_to_flip,'%s_makeup' % sample]

    print("done\n")

    if return_metric:
        return compute_metric(intersection_SNPs, sample = sample, minDP = minDP)
    else:
        return intersection_SNPs



def compute_metric(dat_all_genotypes, sample, minDP = None):
    print('Measure performance...')
    print('Among all variants')

    ### Among all tested, how many are recovered
    N = sum(~dat_all_genotypes['REF_y'].isnull())
    called_percentage = N/float(len(dat_all_genotypes))

    true_hits = np.sum(dat_all_genotypes[sample] == dat_all_genotypes['%s_called' % sample])
    called_correct_percentage = true_hits / float(len(dat_all_genotypes))
    recovered = [len(dat_all_genotypes), N, true_hits, called_percentage, called_correct_percentage]
    print("Among %d variants identified by WGS, %d (%.3f) are called by ATAC-seq reads, %d (%.3f) are correct" % (len(dat_all_genotypes), N, called_percentage, true_hits, called_correct_percentage))


    ### Among the tested variants, evaluate the performance
    print('Among recovered variants by ATAC-seq reads:')
    dat_all_genotypes = dat_all_genotypes[~dat_all_genotypes['REF_y'].isnull()]
    confusion_matrix = obtain_confusion_matrix(dat_all_genotypes, sample)

    recall_arr = np.array(map(float, np.diag(np.array(confusion_matrix))) / np.reshape(np.array(confusion_matrix.sum(axis=1)), [1,3])).ravel()
    precision_arr = np.array(map(float, np.diag(np.array(confusion_matrix))) / np.reshape(np.array(confusion_matrix.sum(axis=0)), [1,3])).ravel()
    performance = list(recall_arr) + list(precision_arr)
    print(["Recall: ", recall_arr])
    print(["Precision: ", precision_arr])
    print("done\n")
    print("================================================\n")

    return [sample, minDP] + recovered + performance




def obtain_confusion_matrix(intersection_df, sample):
    intersection_df = intersection_df.copy()
    intersection_df['Alt1'] = [x.split('/')[0] for x in intersection_df[sample]]
    intersection_df['Alt2'] = [x.split('/')[1] for x in intersection_df[sample]]
    intersection_df['MT']= (intersection_df['Alt1'] != intersection_df['REF_x']) & (intersection_df['Alt2'] != intersection_df['REF_x'])
    intersection_df['HT'] = intersection_df['Alt1'] != intersection_df['Alt2']

    intersection_df['called_Alt1'] = [x.split('/')[0] for x in intersection_df['%s_called' % sample]]
    intersection_df['called_Alt2'] = [x.split('/')[1] for x in intersection_df['%s_called' % sample]]
    intersection_df['called_MT']= (intersection_df['called_Alt1'] != intersection_df['REF_x']) & (intersection_df['called_Alt2'] != intersection_df['REF_x'])
    intersection_df['called_HT'] = intersection_df['called_Alt1'] != intersection_df['called_Alt2']

    ## Real Genotype

    # Heterozygous - AB
    ab = intersection_df[intersection_df['HT']]

    # called genotypes
    ab_called_ab = np.sum(ab['called_HT'])
    ab_called_aa = np.sum((~ab['called_HT']) & (~ab['called_MT']))
    ab_called_bb = np.sum((~ab['called_HT']) & (ab['called_MT']))

    assert (ab_called_ab + ab_called_aa + ab_called_bb == len(ab))

    # Homogenous wildtype
    aa = intersection_df[(~intersection_df['HT']) & (~intersection_df['MT'])]

    # called genotypes
    aa_called_ab = np.sum(aa['called_HT'])
    aa_called_aa = np.sum((~aa['called_HT']) & (~aa['called_MT']))
    aa_called_bb = np.sum((~aa['called_HT']) & (aa['called_MT']))
    assert (aa_called_ab + aa_called_aa + aa_called_bb == len(aa))

    # Homogenous mutation
    bb = intersection_df[(~intersection_df['HT']) & (intersection_df['MT'])]

    # called genotype
    bb_called_ab = np.sum(bb['called_HT'])
    bb_called_aa = np.sum((~bb['called_HT']) & (~bb['called_MT']))
    bb_called_bb = np.sum((~bb['called_HT']) & (bb['called_MT']))
    assert (bb_called_ab + bb_called_aa + bb_called_bb == len(bb))

    confusion_matrix = pd.DataFrame([[aa_called_aa, aa_called_ab, aa_called_bb], 
                                     [ab_called_aa, ab_called_ab, ab_called_bb],
                                     [bb_called_aa, bb_called_ab, bb_called_bb]])

    confusion_matrix.columns = ['call_AA', 'called_AB', 'called_BB']
    confusion_matrix.index = ['true_AA', 'true_AB', 'true_BB']
    
    return confusion_matrix


def readin_INFO(sample):

    sample_info = pd.read_csv('%s/%s.filtered.recode.INFO.vcf' % (variant_calling_dir, sample),sep='\t', header = None)
    sample_info.columns=['#CHROM', 'POS', 'DP', 'PL', 'GQ']
    #sample_info['#CHROM'] = [int(x.replace('chr', '')) for x in sample_info['#CHROM']]

    return sample_info



def compare_DP(sample_df):

    mismatch = sample_df[sample_df[sample] != sample_df['%s_called' % sample]]
    mismatch = mismatch.merge(info_df, on=['#CHROM', 'POS'])

    match = sample_df[sample_df[sample] == sample_df['%s_called' % sample]]
    match = match.merge(info_df, on=['#CHROM', 'POS'])

    x = [np.log10(t) for t in np.array(mismatch['DP'])]
    y = [np.log10(t) for t in np.array(match['DP'])]

    print('Compare minDP for the matched and mismatched SNP sites')
    print([np.mean(np.log10(x))-np.mean(np.log10(y)), ranksums(np.log10(x), np.log10(y))])
    print('\n\n')


