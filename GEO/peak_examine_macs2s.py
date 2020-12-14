import os
import numpy as np
import pandas as pd
import sys
import pdb


def compare_macs2_overlapping(sample1, sample2):
	## Compare peaks from macs2 and Genrich
	samples = pd.read_csv('test_all.txt', header = None)
	samples = np.array(samples[0])

	number_peaks = []
	median = []

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
	print(np.median(median, axis=0))

	number_peaks = pd.DataFrame(number_peaks)
	number_peaks.columns = ['N_peaks_macs2', 'N_peaks_Genrich', 
                        'N_peaks_macs2_in_Genrich', 'N_peaks_Genrich_in_macs2',
                        'bp_macs2', 'bp_Genrich', 'bp_overlap']

	number_peaks.to_csv('evaluate_peaks/Examine_peaks_results_genrich_macs2.txt', sep='\t', index = False)



if __name__ == '__main__':
 	root_dir = sys.argv[1]
        macs_dir = '%s/Peaks' % root_dir
        Genrich_dir = '%s/Peaks_Genrich' % root_dir
        overlap_dir = '%s/Peaks_Examine' % root_dir

	compare_macs2_overlapping(sample1, sample2)





