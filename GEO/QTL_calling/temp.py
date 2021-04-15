import os
import numpy as np
import pandas as pd

import sys
sys.path.append('/work-zfs/abattle4/heyuan/tools/python3_lib/lib/python3.6/site-packages/')

from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
import scipy.stats as stats

import pdb

def compute_RPKM(peaks, bam_stats, peak_dir):
    '''
    Compute RPKM (Reads per kilo base per million mapped reads)
    '''
    print('Compute RPKM')
        
    read_depth = np.array(bam_stats.set_index('Sample').loc[samples]['Reads']) / 1e6
    peak_width = np.array(peaks['END'] - peaks['START']) / 1e3

    x = np.divide(np.array(peaks[samples]) , read_depth)
    RPKM = np.divide(x.transpose(), peak_width).transpose()
    RPKM_dat = pd.DataFrame(RPKM)

    RPKM_dat.columns = samples
    
    RPKM_dat['CHR'] = np.array(peaks['CHR'])
    RPKM_dat['START'] = np.array(peaks['START'])
    RPKM_dat['END'] = np.array(peaks['END'])
    RPKM_dat['PEAK'] = np.array(peaks['PEAK'])
    RPKM_dat = RPKM_dat[['PEAK', 'CHR', 'START', 'END']  + list(samples)]
    RPKM_dat.to_csv('%s/peak_by_sample_matrix_RPKM.txt' % (peak_dir), sep=' ', index = False)
    
    return RPKM_dat


def read_in_peaks(peak_dir, permute):
    print('Read in Peak and convert to RPKM')
    # sample read depth
    try:
        RPKM_dat = pd.read_csv('%s/peak_by_sample_matrix_RPKM.txt' % peak_dir, sep=' ')
    except:
        # Peaks
        peaks = pd.DataFrame()
        for CHROMOSOME in range(1,23):
            peaks_chr = pd.read_csv('%s/peak_by_sample_matrix_chr%d.txt' % (peak_dir, CHROMOSOME), sep=' ')
            peaks_chr = peaks_chr[[x for x in peaks_chr.columns if 'Unnamed' not in x]]
            peaks = peaks.append(peaks_chr)

        RPKM_dat = compute_RPKM(peaks,bam_stats, peak_dir)

    RPKM_dat_counts = RPKM_dat[samples]
    print(RPKM_dat_counts.shape)
 
    if permute:
        values = RPKM_dat_counts.values.ravel()
        permuted_values = np.random.choice(values, replace=False, size = len(values))
        RPKM_dat_counts = pd.DataFrame(permuted_values.reshape(RPKM_dat_counts.shape))
        RPKM_dat_counts.columns = samples

    df = RPKM_dat_counts.copy()
    
    qn_df = normalize_quantiles(df)
    norm_df = inverse_normal_transform(qn_df)
    
    # format
    #samples = [x for x in RPKM_dat.columns if x.startswith('SRR')]

    corrected_dat = pd.DataFrame(norm_df)
    corrected_dat.columns = samples
    corrected_dat['CHR'] = np.array(RPKM_dat['CHR'])
    corrected_dat['START'] = np.array(RPKM_dat['START'])
    corrected_dat['END'] = np.array(RPKM_dat['END'])
    corrected_dat['PEAK'] = np.array(RPKM_dat['PEAK'])

    corrected_dat = corrected_dat[['PEAK', 'CHR', 'START', 'END']  + list(samples)]

    return corrected_dat


# In[4]:


# https://github.com/broadinstitute/gtex-pipeline/blob/b53b734a9b096caed237952b34cbce88b38485bd/qtl/src/rnaseqnorm.py#L53

def normalize_quantiles(df):
    """
    Quantile normalization to the average empirical distribution
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")
    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003
    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    M = df.values.copy()

    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n

    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1

        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1

    return pd.DataFrame(M, index=df.index, columns=df.columns)


def inverse_normal_transform(M):
    """
    Transform rows to a standard normal distribution
    """
    R = stats.mstats.rankdata(M, axis=1)  # ties are averaged
    if isinstance(M, pd.DataFrame):
        Q = pd.DataFrame(stats.norm.ppf(R/(M.shape[1]+1)), index=M.index, columns=M.columns)
    else:
        Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q


def remove_PCs(df, PCA_k = 10):
    print('Fit in PCA')
    df_values = df[samples]
    pca = PCA(n_components=20)
    pca.fit(df_values.transpose()) 
    pcs = pca.transform(df_values.transpose())
    
    reg = LinearRegression().fit(pcs[:, :PCA_k], df_values.transpose())
    dat = df_values - np.dot(reg.coef_, pcs[:, :PCA_k].transpose())

    for col in ['PEAK','CHR', 'START', 'END']:
        dat[col] = df[col]

    return [dat, pca]


if __name__ == '__main__':

    peak_calling = sys.argv[1]

    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/ATAC_seq/alignment_bowtie/'
    PEAK_DIR = '%s/Peaks_%s' % (root_dir, peak_calling)
    BAM_DIR = '%s/first_pass_bqsr' % root_dir

    genotype_available_file = open('%s/VCF_files/merged_vcf_header.vcf' % root_dir, 'r')
    for l in genotype_available_file.readlines():
        if l.startswith('##'):
            continue
        elif l.startswith('#'):
            samples = [x for x in l.rstrip().split('\t') if 'SRR' in x]
        else:
            break

    bam_stats = pd.read_csv('%s/bam_stats.txt' % BAM_DIR, sep=' ', header=None)
    bam_stats.columns=['Sample', 'Reads']
    bam_stats = bam_stats[~bam_stats['Reads'].isnull()]
    bam_stats = bam_stats.drop_duplicates()
    samples = np.intersect1d(samples, bam_stats['Sample'])

    bam_stats = bam_stats.set_index('Sample').loc[samples].reset_index()
 
    centered_data = read_in_peaks(PEAK_DIR, permute = False)
    corrected_data = centered_data.copy()

    for i in range(1,23):
        df_save = corrected_data[corrected_data['CHR'] == i]
        df_save.to_csv('%s/peak_by_sample_matrix_RPKM_corrected_chromosome%d.txt' % (PEAK_DIR, i), sep='\t', index = False)

        df_save['MID'] = list(map(round, (df_save['START'] + df_save['END']) / 2))
        df_save['MID2'] = df_save['MID'] + 1
        df_save['#CHR'] = ['chr%d' % i] * len(df_save)
        df_save = df_save.sort_values('MID')
        df_save[['#CHR', 'MID', 'MID2', 'PEAK'] + list(samples)].to_csv('%s/peak_by_sample_matrix_RPKM_corrected_chromosome%d_forFastQTL.bed'% (PEAK_DIR, i), sep='\t', index=False)

