import os
import numpy as np
import pandas as pd
import sys
import pdb


def compare_macs2_overlapping(sample1, sample2):
	number_peaks = []
	median = []

	macs_peaks_s1 = pd.read_csv('%s/%s_peaks.narrowPeak' % (macs_dir, sample1), sep='\t', header = None)
	macs_peaks_s2 = pd.read_csv('%s/%s_peaks.narrowPeak' % (macs_dir, sample2), sep='\t', header = None)
    
	# look at overlap
	macs2_overlap_peaks = pd.read_csv('%s/%s_%s_macs2_overlap.bed' % (overlap_dir, sample1, sample2), sep='\t', header = None)
	macs2_overlap_peaks = macs2_overlap_peaks.drop_duplicates()
    
	s1_basepair = pd.read_csv('%s/%s_macs2_peak_regions.bed' % (overlap_dir, sample1), sep='\t', header = None)
	s2_basepair = pd.read_csv('%s/%s_macs2_peak_regions.bed' % (overlap_dir, sample2), sep='\t', header = None)
	overlap_bp = pd.read_csv('%s/%s_%s_macs2_overlap.bed' % (overlap_dir, sample1, sample2), sep='\t', header = None)

	number_peaks.append([len(macs_peaks_s1), len(macs_peaks_s2), 
                         len(set(macs2_overlap_peaks[3])), 
                         np.sum(s1_basepair[2] - s1_basepair[1]), 
                         np.sum(s2_basepair[2] - s2_basepair[1]),
                         np.sum(overlap_bp[2] - overlap_bp[1])])

	number_peaks = pd.DataFrame(number_peaks)
	number_peaks.columns = ['N_peaks_s1', 'N_peaks_s2', 
                        'N_peaks_s1_s2',
                        'bp_s1', 'bp_s2', 'bp_overlap']

	number_peaks.to_csv('evaluate_peaks/Examine_peaks_results_macs2_%s_%s.txt' % (sample1, sample2), sep='\t', index = False)



if __name__ == '__main__':
	sample1 = sys.argv[1]
	sample2 = sys.argv[2]
 	root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/PRJNA603385/ATAC_seq/alignment_bowtie/'
        macs_dir = '%s/Peaks' % root_dir
	overlap_dir = '%s/Peaks_Examine' % root_dir

	compare_macs2_overlapping(sample1, sample2)





