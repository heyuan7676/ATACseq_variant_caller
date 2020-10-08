import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ranksums


def read_in_WGS_GT(sample):
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/Genotype'

    chromosome = 22
    WGS_fn = '%s/1k_genome_chr%d.genotypes.tsv' % (WGS_dir, chromosome)
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t', nrows = 10)
    try:
        col_to_read = list(WGS_result.columns).index(sample)
    except:
        print('%s does not have genotype data' % sample)
        return 
    
    WGS_result = pd.DataFrame()

    for chromosome in range(1,23):
        WGS_fn = '%s/1k_genome_chr%d.genotypes.tsv' % (WGS_dir, chromosome)
        token = pd.read_csv(WGS_fn, 
                            comment = '$', 
                            sep='\t', 
                            usecols=[0,1,2,col_to_read])
        WGS_result = WGS_result.append(token)    
    WGS_result = WGS_result.drop_duplicates() 
    return WGS_result



def obtain_intersection_df(sample, minDP, WGS_result, restrict_to_SNP = True):
    
    WGS_result = WGS_result.copy()
    
    SNP_calling_dir = '/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/ATAC_seq/alignment_bowtie'
        
    SNP_called_fn = '%s/%s-RG-dedup-hapcal.filtered.genotype.minDP%d.txt' % (SNP_calling_dir, sample, minDP)
    SNP_called = pd.read_csv(SNP_called_fn, comment="$", sep='\t', nrows = 10)
    try:
        col_to_read = np.where([sample in x for x in SNP_called.columns])[0][0]
    except:
        print('%s does not called genotpye' % sample)
        return
    
    SNP_called = pd.read_csv(SNP_called_fn, comment="$", sep='\t', usecols=[0, 1, 2, col_to_read])
    SNP_called.columns = list(SNP_called.columns[:3]) + ['%s_called' % sample]
    SNP_called['#CHROM'] = [int(x.replace('chr','')) for x in SNP_called['#CHROM']]
    
    # restrict to SNP
    if restrict_to_SNP:
	WGS_result = WGS_result.iloc[np.where([len(x) == 1 for x in WGS_result['REF']])[0]]
        WGS_result = WGS_result.iloc[np.where([len(x.split('/')[0]) == 1 for x in WGS_result[sample]])[0]]
        WGS_result = WGS_result.iloc[np.where([len(x.split('/')[1]) == 1 for x in WGS_result[sample]])[0]]
        
        SNP_called = SNP_called.iloc[np.where([len(x.split('/')[1]) == 1 for x in SNP_called['%s_called' % sample]])[0]]
        SNP_called = SNP_called.iloc[np.where([len(x.split('/')[1]) == 1 for x in SNP_called['%s_called' % sample]])[0]]
    
    # variants with some read
    intersection_SNPs = WGS_result.merge(SNP_called, on=['#CHROM', 'POS'])
    
    # match for the HT order
    intersection_SNPs['%s_makeup' % sample] = ['/'.join(x.split('/')[::-1]) for x in intersection_SNPs[sample]]
    
    HT = np.where(intersection_SNPs[sample] != intersection_SNPs['%s_makeup' % sample])[0]
    match_idx = np.where(intersection_SNPs['%s_called' % sample] == intersection_SNPs['%s_makeup' % sample])[0]
    need_to_flip = np.intersect1d(HT, match_idx)
    intersection_SNPs.loc[need_to_flip,sample] = intersection_SNPs.loc[need_to_flip,'%s_makeup' % sample]

    true_hits = np.sum(intersection_SNPs[sample] == intersection_SNPs['%s_called' % sample])

    print('For %s:' % sample)
    print("In total, WGS called %d SNVs" % len(WGS_result))
    print("Variant caller called %d SNVs" % len(SNP_called))
    print("Variant called called %d (%.3f) SNVs found in WGS" % (len(intersection_SNPs), 
                                                                 len(intersection_SNPs)/len(SNP_called)))
    print("   %d (%.3f) are correct\n" % (true_hits, true_hits / len(SNP_called)))
    
    
    return intersection_SNPs



def obtain_confusion_matrix(intersection_df):
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

    sample_info = pd.read_csv('%s/%s-RG-dedup-hapcal.filtered.variants.recode.INFO.txt' % (variant_calling_dir, sample),sep='\t', header = None)
    sample_info.columns=['#CHROM', 'POS', 'DP', 'PL']
    sample_info['#CHROM'] = [int(x.replace('chr', '')) for x in sample_info['#CHROM']]

    return sample_info



def compare_DP(sample_df):

    mismatch = sample_df[sample_df[sample] != sample_df['%s_called' % sample]]
    mismatch = mismatch.merge(info_df, on=['#CHROM', 'POS'])

    match = sample_df[sample_df[sample] == sample_df['%s_called' % sample]]
    match = match.merge(info_df, on=['#CHROM', 'POS'])

    x = [np.log10(x) for x in np.array(mismatch['DP'])]
    y = [np.log10(x) for x in np.array(match['DP'])]

    print('Compare minDP for the matched and mismatched SNP sites')
    print([np.mean(np.log10(x))-np.mean(np.log10(y)), ranksums(np.log10(x), np.log10(y))])
    print('\n\n')



if __name__ == '__main__':

    sample = sys.argv[1]

    # WGS genotype
    WGS_df = read_in_WGS_GT(sample)

    # variant calling information
    variant_calling_dir = '/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/ATAC_seq/alignment_bowtie'
    info_df = readin_INFO(sample)
    info_df_unamb = info_df[['0,0' not in x for x in info_df['PL']]]

    for minDP in range(3,11):
        print("Use Filter: minDP >= %d" % minDP)
        sample_df = obtain_intersection_df(sample, minDP = minDP, WGS_result=WGS_df)
        if sample_df is not None:
            confusion_matrix = obtain_confusion_matrix(sample_df)
            confusion_matrix.to_csv('confusion_matrices/%s_bowtie_GATK_DPmin%d_SNPs.txt' % (sample, minDP), sep='\t')
            compare_DP(sample_df)
    
            # remove ambiguous calls
            sample_df = sample_df.merge(info_df_unamb, on = ['#CHROM', 'POS'])
            confusion_matrix = obtain_confusion_matrix(sample_df)
            confusion_matrix.to_csv('confusion_matrices/%s_bowtie_GATK_DPmin%d_SNPs_noAmbig.txt' % (sample, minDP), sep='\t')
            
        else:
            continue


