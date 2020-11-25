import os
import numpy as np
import pandas as pd

import sys
sys.path.append('/work-zfs/abattle4/heyuan/tools/python3_lib/lib/python3.6/site-packages/')

from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
import scipy.stats as stats

import pdb

def compute_RPKM(peaks, bam_stats, samples, peak_dir):
    '''
    Compute RPKM (Reads per kilo base per million mapped reads)
    '''
    print('Compute RPKM')
        
    read_depth = np.array(bam_stats.set_index('Sample').loc[peaks.columns[4:]]['Reads']) / 1e6
    peak_width = np.array(peaks['END'] - peaks['START']) / 1e3

    x = np.divide(np.array(peaks[peaks.columns[4:]]) , read_depth)
    RPKM = np.divide(x.transpose(), peak_width).transpose()
    RPKM_dat = pd.DataFrame(RPKM)

    RPKM_dat.columns = samples
    
    RPKM_dat['CHR'] = np.array(peaks['CHR'])
    RPKM_dat['START'] = np.array(peaks['START'])
    RPKM_dat['END'] = np.array(peaks['END'])
    RPKM_dat['PEAK'] = np.array(peaks['PEAK'])
    RPKM_dat = RPKM_dat[['PEAK', 'CHR', 'START', 'END']  + samples]
    RPKM_dat.to_csv('%s/peak_by_sample_matrix_RPKM.txt' % (peak_dir), sep=' ', index = False)
    
    return RPKM_dat


def read_in_peaks(bam_dir, peak_dir, permute):
    bam_stats = pd.read_csv('%s/bam_stats.txt' % bam_dir, sep='\t', header=None)
    bam_stats.columns=['Sample', 'Reads']

    try:
        RPKM_dat = pd.read_csv('%s/peak_by_sample_matrix_RPKM.txt' % (peak_dir), sep=' ')
        samples = [x for x in RPKM_dat.columns if x.startswith('HG')]
        print('Read in RPKM')
                
    except:
        print('Read in Peak and convert to RPKM')
        # Peaks
        peaks = pd.DataFrame()
        for CHROMOSOME in range(1,23):
            peaks_chr = pd.read_csv('%s/peak_by_sample_matrix_chr%d.txt' % (peak_dir, CHROMOSOME), sep=' ')
            peaks_chr = peaks_chr[[x for x in peaks_chr.columns if 'Unnamed' not in x]]
            peaks = peaks.append(peaks_chr)
            
        samples = [x for x in peaks.columns if x.startswith('HG')]

        # sample read depth
        RPKM_dat = compute_RPKM(peaks,bam_stats,samples, peak_dir)
        
    RPKM_dat_counts = RPKM_dat[samples]
    
    if permute:
        values = RPKM_dat_counts.values.ravel()
        permuted_values = np.random.choice(values, replace=False, size = len(values))
        RPKM_dat_counts = pd.DataFrame(permuted_values.reshape(RPKM_dat_counts.shape))
        RPKM_dat_counts.columns = samples

    df = RPKM_dat_counts.copy()
    
    qn_df = normalize_quantiles(df)
    norm_df = inverse_normal_transform(qn_df)
    
    # format
    samples = [x for x in RPKM_dat.columns if x.startswith('HG')]

    corrected_dat = pd.DataFrame(norm_df)
    corrected_dat.columns = samples
    corrected_dat['CHR'] = np.array(RPKM_dat['CHR'])
    corrected_dat['START'] = np.array(RPKM_dat['START'])
    corrected_dat['END'] = np.array(RPKM_dat['END'])
    corrected_dat['PEAK'] = np.array(RPKM_dat['PEAK'])

    corrected_dat = corrected_dat[['PEAK', 'CHR', 'START', 'END']  + samples]

    return [corrected_dat, bam_stats, samples]


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


def remove_PCs(df, samples, PCA_k = 10):

    df_values = df[samples]
    pca = PCA(n_components=50)
    pca.fit(df_values.transpose()) 
    pcs = pca.transform(df_values.transpose())
    
    reg = LinearRegression().fit(pcs[:, :PCA_k], df_values.transpose())
    dat = df_values - np.dot(reg.coef_, pcs[:, :PCA_k].transpose())

    for col in ['PEAK','CHR', 'START', 'END']:
        dat[col] = df[col]

    return [dat, pca]


if __name__ == '__main__':

    #PEAK_DIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_Genrich/combined/'
    #BAM_DIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/first_pass_bqsr'

    PEAK_DIR = sys.argv[1]
    BAM_DIR = sys.argv[2]

    [centered_data, bam_DF, samples] = read_in_peaks(BAM_DIR, PEAK_DIR, permute = False)
    [corrected_data, pca_model] = remove_PCs(centered_data, samples, PCA_k=12)

    [permuted_centered_data, bam_DF, samples] = read_in_peaks(BAM_DIR, PEAK_DIR, True)
    [permuted_corrected_data, permuted_pca_model] = remove_PCs(permuted_centered_data, samples, PCA_k=10)

    idx = np.where( pca_model.explained_variance_ratio_[:50] <  permuted_pca_model.explained_variance_ratio_[:50])[0][0]
    print('Correct for %d PCs' % idx)

    [corrected_data, pca_model] = remove_PCs(centered_data, samples, PCA_k=idx - 1)
    for i in range(1,23):
        df_save = corrected_data[corrected_data['CHR'] == i]
        df_save.to_csv('%s/peak_by_sample_matrix_RPKM_corrected_chromosome%d.txt' % (PEAK_DIR, i), sep='\t', index = False)
        #df_save[['PEAK'] + samples].to_csv('%s/peak_by_sample_matrix_RPKM_corrected_chromosome%d_values.txt' % (PEAK_DIR, i), sep='\t', index = False)
        #df_save[['PEAK','CHR', 'START', 'END']].to_csv('%s/peak_by_sample_matrix_RPKM_corrected_chromosome%d_loc.bed'% (PEAK_DIR, i), sep='\t', index=False)

    #pca_save = pcs[:, :PCA_k].transpose()
    #pca_save = pd.DataFrame(pca_save)
    #pca_save.index = ['PC%d' % i for i in range(1, PCA_k+1)]
    #pca_save.columns = samples
    #pca_save.to_csv('%s/PCs.txt' % PEAK_DIR, sep='\t')

    


