import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys
import pdb



def compute_TPM(peaks, bam_stats, samples, peak_dir):
    '''
    Compute TPM (Transcripts Per Kilobase Million)
    '''
    print('Compute TPM')

    read_depth = np.array(bam_stats.set_index('Sample').loc[peaks.columns[4:]]['Reads']) / 1e6
    peak_width = np.array(peaks['END'] - peaks['START']) / 1e3

    x = np.divide(np.array(peaks[peaks.columns[4:]]).transpose(), peak_width).transpose()
    TPM = np.divide(x, read_depth)
    TPM_dat = pd.DataFrame(TPM)

    TPM_dat.columns = samples

    TPM_dat['CHR'] = np.array(peaks['CHR'])
    TPM_dat['START'] = np.array(peaks['START'])
    TPM_dat['END'] = np.array(peaks['END'])
    TPM_dat['PEAK'] = np.array(peaks['PEAK'])
    TPM_dat = TPM_dat[['PEAK', 'CHR', 'START', 'END']  + samples]
    TPM_dat.to_csv('%s/peak_by_sample_matrix_TPM.txt' % (peak_dir), sep=' ', index = False)

    return TPM_dat



def read_in_peaks_TPM(bam_dir, peak_dir, permute):
    bam_stats = pd.read_csv('%s/bam_stats.txt' % bam_dir, sep='\t', header=None)
    bam_stats.columns=['Sample', 'Reads']

    try:
        TPM_dat = pd.read_csv('%s/peak_by_sample_matrix_TPM.txt' % (peak_dir), sep=' ')
        samples = [x for x in TPM_dat.columns if x.startswith('HG')]
        print('Read in TPM')
                
    except:
        print('Read in Peak and convert to TPM')
        # Peaks
        peaks = pd.DataFrame()
        for CHROMOSOME in range(1,23):
            peaks_chr = pd.read_csv('%s/peak_by_sample_matrix_chr%d.txt' % (peak_dir, CHROMOSOME), sep=' ')
            peaks_chr = peaks_chr[[x for x in peaks_chr.columns if 'Unnamed' not in x]]
            peaks = peaks.append(peaks_chr)
           
        samples = [x for x in peaks.columns if x.startswith('HG')]

        # sample read depth
        TPM_dat = compute_TPM(peaks,bam_stats,samples, peak_dir)
        


if __name__ == '__main__':

    peak_calling = sys.argv[1]

    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/'

    if peak_calling == 'macs2':
        PEAK_DIR = '%s/Peaks_MACS2' % root_dir

    if peak_calling == 'macs2_combined':
        PEAK_DIR = '%s/Peak_version1/Peaks/combined' % root_dir

    elif peak_calling == 'Genrich':
        PEAK_DIR = '%s/Peaks_Genrich' % root_dir

    elif peak_calling == 'Genrich_combined':
        PEAK_DIR = '%s/Peaks_Genrich/combined' % root_dir

    BAM_DIR = '%s/first_pass_bqsr' % root_dir

    [centered_data, bam_DF, samples] = read_in_peaks_TPM(BAM_DIR, PEAK_DIR, permute = False)
