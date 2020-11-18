import os
import numpy as np
import pandas as pd
import sys

## Compare peaks from macs2 and Genrich

macs_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks'
Genrich_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_Genrich'
overlap_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_Examine'


## Compare peaks from macs2 and Genrich
samples = pd.read_csv('test_all.txt', header = None)
samples = np.array(samples[0])

median = []
number_peaks = []

for s in samples:
    macs_peaks = pd.read_csv('%s/%s_peaks.narrowPeak' % (macs_dir, s), sep='\t', header = None)
    Genrich_peaks = pd.read_csv('%s/%s_peaks.Genrich.narrowPeak' % (Genrich_dir, s), sep='\t', header = None)
    median.append([np.median(np.array(macs_peaks[2] - macs_peaks[1])),
                   np.median(np.array(Genrich_peaks[2] - Genrich_peaks[1]))])
    
    # look at overlap
    macs2_in_genrich = pd.read_csv('%s/%s_macs2_in_genrich.bed' % (overlap_dir, s), sep='\t', header = None)
    macs2_in_genrich = macs2_in_genrich.drop_duplicates()
    genrich_in_macs2 = pd.read_csv('%s/%s_genrich_in_macs2.bed' % (overlap_dir, s), sep='\t', header = None)
    genrich_in_macs2 = genrich_in_macs2.drop_duplicates()
    
    # overlapping basepairs
    macs2_basepair = pd.read_csv('%s/%s_macs2_peak_regions.bed' % (overlap_dir, s), sep='\t', header = None)
    genrich_basepair = pd.read_csv('%s/%s_genrich_peak_regions.bed' % (overlap_dir, s), sep='\t', header = None)
    overlap_bp = pd.read_csv('%s/%s_macs2_genrich_overlap.bed' % (overlap_dir, s), sep='\t', header = None)
    number_peaks.append([len(macs_peaks), len(Genrich_peaks), 
                         len(set(macs2_in_genrich[3])), 
                         len(set(genrich_in_macs2[3])),
                         np.sum(macs2_basepair[2] - macs2_basepair[1]), 
                         np.sum(genrich_basepair[2] - genrich_basepair[1]),
                         np.sum(overlap_bp[2] - overlap_bp[1])])
    
median = np.array(median)
number_peaks = pd.DataFrame(number_peaks)
number_peaks.columns = ['N_peaks_macs2', 'N_peaks_Genrich', 
                        'N_peaks_macs2_in_Genrich', 'N_peaks_Genrich_in_macs2',
                        'bp_macs2', 'bp_Genrich', 'bp_overlap']

number_peaks.to_csv('Examine_peaks_results.txt', sep='\t', index = False)


