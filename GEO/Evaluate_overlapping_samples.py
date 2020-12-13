import pdb
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ranksums
sys.path.append('/work-zfs/abattle4/heyuan/Variant_calling/scripts/QTL')
sys.path.append('/work-zfs/abattle4/heyuan/Variant_calling/GBR_benchmarking/Genotype_calling/Evaluation')
from prepare_data_matrix import achieve_ll
#from Evaluation_metrics import *


def read_in_1k_variants(chromosome):
    fn = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR//Genotype/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.variants.txt' % chromosome
    wgs_1k = pd.read_csv(fn, usecols = [1, 2], sep=' ')
    return wgs_1k



def obtain_GT_df(sample, oneK_variants, restrict_to_SNP = True, return_df = True, Imputed = False, minDP = 2):

    print('Read in genotype data called from ATAC-seq reads...')

    if Imputed:
        SNP_calling_dir = '%s/Imputation/minDP%d/%s' % (root_dir, minDP, sample)
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
    SNP_called = SNP_called[~SNP_called[['#CHROM', 'POS']].duplicated()]

    SNP_called['#CHROM'] = SNP_called['#CHROM'].apply(str)
    oneK_variants['#CHROM'] = oneK_variants['#CHROM'].apply(str)

    SNP_called_in1k = SNP_called.merge(oneK_variants, on = ["#CHROM", "POS"])
    print('Keep variants that have MAF >= 0.05 from 1000 Genome project: %d --> %d variants' % (len(SNP_called), len(SNP_called_in1k)))
    SNP_called = SNP_called_in1k.copy()


    # restrict to SNP
    if restrict_to_SNP:
        SNP_called = SNP_called.iloc[np.where([len(x.split('/')[1]) == 1 for x in SNP_called['%s_called' % sample]])[0]]
        SNP_called = SNP_called.iloc[np.where([len(x.split('/')[1]) == 1 for x in SNP_called['%s_called' % sample]])[0]]


    return SNP_called


def convert_gt_to_number(arri):
    '''
    convert genotype type A/A to 0, 1, 2 ...
    '''
    arri_idx = np.where(arri!='0')[0]
    nts = np.unique([a for b in [x.split('/') for x in arri[arri_idx]] for a in b])
    # remove the sites with more than two alleles
    if len(nts) > 2:
        return np.ones(len(arri)) * (-1)
    nts_dict = {}
    i = 0
    for ni in nts:
        nts_dict[ni] = i
        i += 1
    numeric_gt = [np.sum([nts_dict[a] for a in x.split('/')]) for x in arri[arri_idx]]
    numeric_gt_arr = np.ones(len(arri)) * (-1)
    numeric_gt_arr[arri_idx] = numeric_gt
    return numeric_gt_arr



def readin_genotype_called_and_imputed(sample, include_imputation = False, restrict_to_SNP = False):
    if restrict_to_SNP:
        suffix = '_allSNPs'
    else:
	suffix = '_allVariants'

    # read in variants from 1k Genome
    oneK_snps = pd.DataFrame()
    for c in range(1, 23):
	snps = read_in_1k_variants(c)
	oneK_snps = oneK_snps.append(snps)

    # variant calling information
    print('Evaluate genotype originally called from ATAC-seq reads...')
    orginally_called = obtain_GT_df(sample, oneK_snps, restrict_to_SNP=restrict_to_SNP, Imputed = False, minDP = minDP)
    weight_fn = '%s/VCF_files/%s.filtered.recode.INFO.formatted.vcf' % (root_dir, sample)
    try:
        weight_dat = pd.read_csv('%s_PPnumbers' % weight_fn, sep='\t')
    except:
        weight_dat = pd.read_csv(weight_fn, sep=' ', header = None)
        weight_dat.columns=['CHR_POS', 'PL']
        weight_dat = weight_dat.copy()
        weight_dat['PP'] = achieve_ll(np.array(weight_dat['PL']))
        weight_dat.to_csv('%s_PPnumbers' % weight_fn, sep='\t', index = False)

    if include_imputation:
        print('Evaluate genotype imputed')
        imputed = obtain_GT_df(sample, oneK_snps, restrict_to_SNP=restrict_to_SNP, Imputed = True, minDP = minDP)

        GT_DAT = orginally_called.merge(imputed, on = ['#CHROM', 'POS'], how = 'outer')
        GT_DAT = GT_DAT[~GT_DAT['REF_x_x'].isnull()]
    else:
	GT_DAT = orginally_called.copy()

    return GT_DAT





if __name__ == '__main__':
    minDP = int(sys.argv[1])
    sample1 = sys.argv[2]
    sample2 = sys.argv[3]
    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/PRJNA603385/ATAC_seq/alignment_bowtie'
    if not os.path.exists('performance/combined'):
        os.makedirs('performance/combined')

    df1 = readin_genotype_called_and_imputed(sample1, restrict_to_SNP = False)
    df2 = readin_genotype_called_and_imputed(sample2, restrict_to_SNP = False)

    df1_df2 = df1.merge(df2, on = ['#CHROM', 'POS'], how = 'outer')
    only_df1 = np.sum(~df1_df2['%s_called' % sample1].isnull())
    only_df2 = np.sum(~df1_df2['%s_called' % sample2].isnull())

    df1_df2 = df1_df2[~df1_df2['%s_called' % sample1].isnull()]
    df1_df2 = df1_df2[~df1_df2['%s_called' % sample2].isnull()]
    true_gt = map(convert_gt_to_number, np.array(df1_df2[['%s_called' % sample1, '%s_called' % sample2]]))

    true_gt = np.array(true_gt)
    consistent = np.sum(true_gt[:, 0] == true_gt[:, 1])

    inconsistent = true_gt[true_gt[:, 0] != true_gt[:, 1]]
    HT_conversion = np.sum([(x[0] == 1) | (x[1] == 1) for x in inconsistent])

    metrics = [sample1, sample2, only_df1, only_df2, len(df1_df2), consistent, HT_conversion, len(inconsistent) - HT_conversion]
    print(metrics)
	













