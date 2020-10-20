import numpy as np
import pandas as pd
import os

import sys

def read_in_WGS_GT():

    WGS_fn = '%s/1k_genome_chr%d.genotypes.tsv' % (WGS_dir, CHROMOSOME)
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t')
    WGS_result = WGS_result.drop_duplicates()
    
    samples = [x for x in WGS_result.columns if x.startswith('HG')]
    WGS_result.index = WGS_result[['#CHROM', 'POS']].apply(lambda x: '_'.join((str(x[0]), str(x[1]))), axis=1)
    
    # use only the SNPs from called genotypes
    QTL_by_called_genotypes = pd.read_csv('%s/CHR%d_acQTLs_WINDOW_%dkb.txt' % (QTL_dir, CHROMOSOME, WINDOW/1000), sep='\t')    
    WGS_result = WGS_result.loc[QTL_by_called_genotypes['CHR_POS']]
    
    WGS_result = WGS_result.drop_duplicates()
    WGS_result.columns = ['CHR', 'POS'] + list(WGS_result.columns[2:])
    WGS_result.index.name = 'CHR_POS'
    WGS_result = WGS_result.reset_index()
    
    return WGS_result



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
    nts_dict = {}
    i = 0
    for ni in nts:
        nts_dict[ni] = i
        i += 1
    numeric_gt = [np.sum([nts_dict[a] for a in x.split('/')]) for x in arri[arri_idx]]
    numeric_gt_arr = np.ones(len(arri)) * (-1)
    numeric_gt_arr[arri_idx] = numeric_gt
    return numeric_gt_arr




def obtain_numerical_gt(gt_dat):
    '''
    Very slow
    '''
    samples = [x for x in gt_dat.columns if x.startswith('HG')]
    
    gt_numerical_dat = []
    for i in range(len(gt_dat)):
        arri = np.array(gt_dat[samples].iloc[i])
        token = convert_gt_to_number(arri)
        gt_numerical_dat.append(token)
        if i % 500 == 0:
            print('%d rows finished' % i)
            
    gt_numerical_dat = pd.DataFrame(gt_numerical_dat)
    gt_numerical_dat.columns = samples
    gt_numerical_dat['CHR_POS'] = gt_dat['CHR_POS']
    gt_numerical_dat['CHR'] = gt_dat['CHR']
    gt_numerical_dat['POS'] = gt_dat['POS']
    gt_numerical_dat = gt_numerical_dat[['CHR_POS', 'CHR', 'POS']  + samples]
    
    gt_numerical_dat.to_csv('%s/gt_by_sample_matrix_chr%d_numericalGT_Window%dkb_Real.txt' % (Genotype_dir, WINDOW/1000, CHROMOSOME), sep=' ', index = False)
    
    return gt_numerical_dat
    



if __name__ == "__main__":
    CHROMOSOME = int(sys.argv[1])
    WINDOW = int(sys.argv[2])

    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'
    PEAK_DIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks'
    BAM_DIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/first_pass_bqsr'
    Genotype_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Called_GT/'
    VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files'
    QTL_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/QTLs'

    [PEAK_DAT, samples] = read_in_peaks()
    try:
        gt_numerical_dat = pd.read_csv('%s/gt_by_sample_matrix_chr%d_numericalGT_Window%dkb_Real.txt' % (Genotype_dir, WINDOW/1000, CHROMOSOME), sep=' ')
    except:
        print('Convert genotype to numerical values')
        WGS_dat = read_in_WGS_GT()
        WGS_dat = WGS_dat[list(WGS_dat.columns[:3]) + list(np.intersect1d(samples, WGS_dat.columns))]
        WGS_DAT = obtain_numerical_gt(WGS_dat)


