import os
import numpy as np
import pandas as pd

import sys
sys.path.append('/work-zfs/abattle4/heyuan/tools/python3_lib/lib/python3.6/site-packages/')

from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import scipy.stats as stats

import pdb



def scale_expressed_peaks():
    try:
        TPM_dat = pd.read_csv('%s/peak_by_sample_matrix_TPM_expressed_Peaks.txt' % (PEAK_DIR), sep = '\t')
        return TPM_dat
    except:
        TPM_dat = pd.DataFrame()
        for chromosome in range(1,23):
            tpm_dat_chr = pd.read_csv('%s/peak_by_sample_matrix_TPM_chr%d.txt' % (PEAK_DIR, chromosome), sep=' ')
            samples = [x for x in tpm_dat_chr.columns if x.startswith('SRR')]
	    samples = list(set(samples) - set(['SRR10075846', 'SRR5128067', 'SRR8171284']))
            # remove peaks with low TPM
            samples_with_high_tpm = np.sum(tpm_dat_chr[samples]> 0.4, axis=1)
            high_tpm = np.where(samples_with_high_tpm > len(samples) / 5)[0]
            print(chromosome, len(high_tpm))
            tpm_dat_chr = tpm_dat_chr.iloc[high_tpm]
            TPM_dat = TPM_dat.append(tpm_dat_chr)

        print('Scale the features')
        TPM_dat = TPM_dat.reset_index(drop = True)
        df_values = TPM_dat[samples].transpose()  # sample x feature
        df_values = scale(df_values, axis = 0)  # scale each feature
        TPM_dat[samples] = df_values.transpose()

        TPM_dat.to_csv('%s/peak_by_sample_matrix_TPM_expressed_Peaks.txt' % (PEAK_DIR), sep='\t', index = False)

        return TPM_dat




def derive_PCs():
    try:
        pcs = pd.read_csv('%s/peak_by_sample_matrix_TPM_topPeak_PCs.txt' % (PEAK_DIR), sep='\t')
        return pcs
    except:
        print('Derive PCs')

    TPM_dat = scale_expressed_peaks()
    samples = [x for x in TPM_dat.columns if x.startswith('SRR')]
    TPM_dat_counts = TPM_dat[samples]

    # select top 500 peaks with highest variance
    peak_sd =  TPM_dat_counts.apply(lambda x: np.std(x), axis=1)
    top_peaks = np.argsort(peak_sd)[::-1][:1000]
    TPM_dat_counts_for_PCA = TPM_dat_counts.iloc[top_peaks]

    # compute PCA
    print('Fit in PCA')

    pca = PCA(n_components=50)
    pca.fit(df_values)
    pcs = pca.transform(df_values)
    pcs = pd.DataFrame(pcs)
    pcs.index = samples
    pcs.to_csv('%s/peak_by_sample_matrix_TPM_topPeak_PCs.txt' % (PEAK_DIR), sep='\t')

    print(pca.explained_variance_ratio_)

    return pcs



def correct_for_PCs(k = 5):
    TPM_dat = scale_expressed_peaks()
    pcs = derive_PCs()
    samples = [x for x in TPM_dat.columns if x.startswith('SRR')]
    df_values = TPM_dat[samples]

    # remove PCs
    reg = LinearRegression().fit(pcs[:, :PCA_k], df_values.transpose())
    corrected_dat = df_values - np.dot(reg.coef_, pcs[:, :PCA_k].transpose())

    for col in ['PEAK','CHR', 'START', 'END']:
        corrected_dat[col] = TPM_dat[col]

    corrected_dat = corrected_dat[['PEAK', 'CHR', 'START', 'END']  + samples]

    return corrected_dat




if __name__ == '__main__':

    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove'
    PEAK_DIR = '%s/Peaks_MACS2' % root_dir

    TPM = scale_expressed_peaks()
    PCs = derive_PCs()






