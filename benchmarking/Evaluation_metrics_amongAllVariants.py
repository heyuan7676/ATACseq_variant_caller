import pdb
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ranksums


def read_in_WGS_GT(sample):
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'

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



def obtain_atac_variants_df(sample, minDP, WGS_result, restrict_to_SNP = True):
    
    WGS_result = WGS_result.copy()
    
    SNP_calling_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Called_GT'
        
    SNP_called_fn = '%s/%s.filtered.genotype.minDP%d.txt' % (SNP_calling_dir, sample, minDP)
    SNP_called = pd.read_csv(SNP_called_fn, comment="$", sep='\t', nrows = 10)
    try:
        col_to_read = np.where([sample in x for x in SNP_called.columns])[0][0]
    except:
        print('%s does not called genotpye' % sample)
        return
    
    SNP_called = pd.read_csv(SNP_called_fn, comment="$", sep='\t', usecols=[0, 1, 2, col_to_read])
    SNP_called.columns = list(SNP_called.columns[:3]) + ['%s_called' % sample]
    #SNP_called['#CHROM'] = [int(x.replace('chr','')) for x in SNP_called['#CHROM']]
    
    # restrict to SNP
    if restrict_to_SNP:
	WGS_result = WGS_result.iloc[np.where([len(x) == 1 for x in WGS_result['REF']])[0]]
        WGS_result = WGS_result.iloc[np.where([len(x.split('/')[0]) == 1 for x in WGS_result[sample]])[0]]
        WGS_result = WGS_result.iloc[np.where([len(x.split('/')[1]) == 1 for x in WGS_result[sample]])[0]]
        
        SNP_called = SNP_called.iloc[np.where([len(x.split('/')[1]) == 1 for x in SNP_called['%s_called' % sample]])[0]]
        SNP_called = SNP_called.iloc[np.where([len(x.split('/')[1]) == 1 for x in SNP_called['%s_called' % sample]])[0]]
   
    # variants with some read
    intersection_SNPs = WGS_result.merge(SNP_called, on=['#CHROM', 'POS'], how = 'left')
    assert len(intersection_SNPs) == len(WGS_result)
    intersection_SNPs.loc[intersection_SNPs['%s_called' % sample].isnull(),'%s_called' % sample] = 'N/N'
    
    # match for the HT order
    intersection_SNPs['%s_makeup' % sample] = ['/'.join(x.split('/')[::-1]) for x in intersection_SNPs[sample]]
    
    HT = np.where(intersection_SNPs[sample] != intersection_SNPs['%s_makeup' % sample])[0]
    match_idx = np.where(intersection_SNPs['%s_called' % sample] == intersection_SNPs['%s_makeup' % sample])[0]
    need_to_flip = np.intersect1d(HT, match_idx)
    intersection_SNPs.loc[need_to_flip,sample] = intersection_SNPs.loc[need_to_flip,'%s_makeup' % sample]

    print('For %s:' % sample)
    N = len(intersection_SNPs) - sum(intersection_SNPs['REF_y'].isnull())
    called_percentage = N/float(len(WGS_result))
    print("Among %d variants identified by WGS, %d (%.3f) are called by ATAC-seq reads" % (len(WGS_result), N, N/float(len(WGS_result))))

    true_hits = np.sum(intersection_SNPs[sample] == intersection_SNPs['%s_called' % sample])
    called_correct_percentage = true_hits / float(len(WGS_result))
    print("   %d (%.3f) are correct\n" % (true_hits, true_hits / float(len(WGS_result))))
    
    return [len(WGS_result), N, true_hits, called_percentage, called_correct_percentage]


def retrive_coverage_for_all_samples(minDP = 2, restrict_to_SNP=True, saveFile='confusion_matrices/Recovered_percentage_SNPs.txt'):
    samples = pd.read_csv('samples.txt', header=None)
    samples = np.array(samples[0])

    recovered = []
    for s in samples[:2]:
        # WGS genotype
        WGS_df = read_in_WGS_GT(s)

        # variant calling information
        print("Use Filter: minDP >= %d" % minDP)
        try:
            [n1, n2, n3, c1,c2] = obtain_atac_variants_df(s, minDP = minDP, WGS_result=WGS_df, restrict_to_SNP=restrict_to_SNP)
            recovered.append([s, n1, n2, n3, c1, c2])
        except:
            print('%s does not have genotype data' % s)

    recovered = pd.DataFrame(recovered)
    recovered.columns = ['Sample', 'N_variants_by_WGS', 'N_variants_by_ATAC', 'N_overlap_variants_by_ATAC', 'Recovered_percentage', 'Correctly_recovered_percentage']
    recovered.to_csv(saveFile, sep='\t', index = False)

    return recovered




if __name__ == '__main__':
    retrive_coverage_for_all_samples()
    #retrive_coverage_for_all_samples(restrict_to_SNP=False, saveFile='confusion_matrices/Recovered_percentage_allVariants.txt')
