import sys
import time
import os
import pdb
import numpy as np
import pandas as pd


def read_in_peaks(PEAK_dir):
    peak_dat = pd.read_csv('%s/peak_by_sample_matrix_RPKM_corrected.txt' % (PEAK_dir), sep='\t', nrows=10)
    samples = [x for x in peak_dat.columns if x.startswith('HG')]
    return samples


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
    


def readin_genotype(Genotype_dir, chromosome, samples, suffix):
    ## Read in Genotpye
    gt_dat = pd.read_csv('%s/gt_by_sample_matrix_chr%d%s.txt' % (Genotype_dir, chromosome, suffix), sep=' ')
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




def readin_genotype_info(gt_dat, VCF_dir, chromosome, samples, suffix):
    ## Read in Genotpye INFO
    info_dat = pd.read_csv('%s/gt_info_by_sample_matrix_chr%d%s.txt' % (VCF_dir, chromosome, suffix), sep=' ')
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


def read_in_WGS_GT(samples, snps = None):
    WGS_fn = '%s/1k_genome_chr%d.genotypes.tsv' % (WGS_dir, CHROMOSOME)
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t')
    WGS_result = WGS_result.drop_duplicates()

    samples = [x for x in WGS_result.columns if x.startswith('HG')]
    WGS_result.index = WGS_result[['#CHROM', 'POS']].apply(lambda x: '_'.join((str(x[0]), str(x[1]))), axis=1)

    if snps is not None:
        # use only the SNPs from called genotypes
        WGS_result = WGS_result.loc[snps]

    WGS_result = WGS_result.drop_duplicates()
    WGS_result.columns = ['CHR', 'POS'] + list(WGS_result.columns[2:])
    WGS_result.index.name = 'CHR_POS'
    WGS_result = WGS_result.reset_index()

    samples = list(np.intersect1d(samples, WGS_result.columns))
    WGS_result = WGS_result[list(WGS_result.columns[:3]) + samples]

    return WGS_result



if __name__ == "__main__":
    CHROMOSOME = int(sys.argv[1])
    #SUFFIX = sys.argv[2]
    GT_DIR = 'minDP2'
    SUFFIX = '_minDP2'
    #SUFFIX = '_GQ1'

    PEAK_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks'
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'

    # root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_subsample_0.5'
    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie'
    Genotype_dir = '%s/Called_GT/%s' % (root_dir, GT_DIR)
    VCF_dir = '%s/VCF_files' % root_dir
    QTL_dir = '%s/QTLs' % root_dir
    #os.makedirs(QTL_dir, exist_ok = True)

    ## read in data
    SAMPLES = read_in_peaks(PEAK_dir)
    GT_DAT = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix=SUFFIX)
    WEIGHT_DAT = readin_genotype_info(gt_dat=GT_DAT, VCF_dir = VCF_dir, chromosome=CHROMOSOME, samples=SAMPLES, suffix=SUFFIX)

    WGS_dat = read_in_WGS_GT(SAMPLES, GT_DAT['CHR_POS'])
    SAMPLES = list(np.intersect1d(SAMPLES, WGS_dat.columns))
    WGS_DAT = obtain_numerical_gt(WGS_dat, SAMPLES)

    WGS_dat_all = read_in_WGS_GT(SAMPLES)
    WGS_DAT_all = obtain_numerical_gt(WGS_dat_all, SAMPLES)


    save_matrix = False
    if save_matrix:
        save_dir = os.path.join(Genotype_dir, 'Save_Matrix')
	os.makedirs(save_dir, exist_ok = True)
    	GT_DAT[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d%s.txt' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index = False)
    	GT_DAT[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d%s_loc.bed' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index=False)

    	WEIGHT_DAT[SAMPLES].to_csv('%s/called_genotypes_chromosome%d%s_weights.txt' % (save_dir, CHROMOSOME, SUFFIX), sep='\t')

    	WGS_DAT[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d%s_realGT.txt' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index = False)
    	WGS_DAT[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d%s_loc_realGT.bed' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index=False)

    	WGS_DAT_all[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d%s_realGT_all.txt' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index = False)
    	WGS_DAT_all[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d%s_loc_realGT_all.bed' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index=False)

    

