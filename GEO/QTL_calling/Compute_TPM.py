import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys
import pdb



def compute_TPM(peaks, peak_dir):
    '''
    Compute TPM (Transcripts Per Kilobase Million)
    '''
    print('Compute TPM')

    bam_stats = pd.read_csv('%s/bam_stats.txt' % BAM_DIR, sep=' ', header=None)
    bam_stats.columns=['Sample', 'Reads']
    bam_stats = bam_stats[~np.isnan(bam_stats['Reads'])]

    samples = [x for x in peaks.columns if x.startswith('SRR')]
    samples = np.intersect1d(np.array(samples), bam_stats['Sample'])

    bam_stats = bam_stats.set_index('Sample').loc[samples]
    peaks = peaks[list(peaks.columns[:4]) + list(samples)]

    read_depth = np.array(bam_stats['Reads']) / 1e6
    peak_width = np.array(peaks['END'] - peaks['START']) / 1e3
   
    x = np.divide(np.array(peaks[peaks.columns[4:]]).transpose(), peak_width).transpose()
    TPM = np.divide(x, read_depth)
    TPM_dat = pd.DataFrame(TPM)

    TPM_dat.columns = samples

    TPM_dat['CHR'] = np.array(peaks['CHR'])
    TPM_dat['START'] = np.array(peaks['START'])
    TPM_dat['END'] = np.array(peaks['END'])
    TPM_dat['PEAK'] = np.array(peaks['PEAKID'])
    TPM_dat = TPM_dat[['PEAK', 'CHR', 'START', 'END']  + list(samples)]
    TPM_dat.to_csv('%s/peak_by_sample_matrix_TPM_chr%d.txt' % (peak_dir, chromosome), sep=' ', index = False)

    return TPM_dat



def read_in_peaks_TPM(bam_dir, peak_dir, permute):
    if 1:
        print('Read in Peak and convert to TPM')
        # Peaks
        peaks = pd.read_csv('%s/overall_matrix_chr%d.txt' % (peak_dir, chromosome), sep='\t')

        # sample read depth
        TPM_dat = compute_TPM(peaks,peak_dir)

    return TPM_dat
        


if __name__ == '__main__':

    chromosome = int(sys.argv[1])

    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove'
    PEAK_DIR = '%s/Peaks_MACS2' % root_dir
    BAM_DIR = '%s/bam_files' % root_dir
    TPM_dat = read_in_peaks_TPM(BAM_DIR, PEAK_DIR, permute = False)
