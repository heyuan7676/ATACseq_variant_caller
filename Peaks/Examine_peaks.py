import os
import numpy as np
import pandas as pd
import sys
import pdb


def compare_macs2_genrich():
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


def compare_macs2_genrich_union():
        median = []
        number_peaks = []

        macs_peaks = pd.read_csv('%s/union-peaks.bed' % (macs_dir), sep='\t', header = None)
        genrich_peaks = pd.read_csv('%s/union-peaks.bed' % (Genrich_dir), sep='\t', header = None)
        median.append([np.median(np.array(macs_peaks[2] - macs_peaks[1])), np.median(np.array(genrich_peaks[2] - genrich_peaks[1]))])

        # look at overlap
        macs2_in_genrich = pd.read_csv('%s/macs2_in_genrich_union.bed' % (overlap_dir), sep='\t', header = None)
        macs2_in_genrich = macs2_in_genrich.drop_duplicates()
        genrich_in_macs2 = pd.read_csv('%s/genrich_in_macs2_union.bed' % (overlap_dir), sep='\t', header = None)
        genrich_in_macs2 = genrich_in_macs2.drop_duplicates()

        # overlapping basepairs
        macs2_basepair = pd.read_csv('%s/macs2_peak_regions_union.bed' % (overlap_dir), sep='\t', header = None)
        genrich_basepair = pd.read_csv('%s/genrich_peak_regions_union.bed' % (overlap_dir), sep='\t', header = None)
        overlap_bp = pd.read_csv('%s/macs2_genrich_overlap_union.bed' % (overlap_dir), sep='\t', header = None)

        number_peaks.append([len(macs_peaks), len(genrich_peaks), len(set(macs2_in_genrich[3])), len(set(genrich_in_macs2[3])), np.sum(macs2_basepair[2] - macs2_basepair[1]), np.sum(genrich_basepair[2] - genrich_basepair[1]), np.sum(overlap_bp[2] - overlap_bp[1])])
        median = np.array(median)
        print(median)

        number_peaks = pd.DataFrame(number_peaks)
        number_peaks.columns = ['N_peaks_macs2', 'N_peaks_genrich', 'N_peaks_macs2_in_genrich', 'N_peaks_genrich_in_macs2', 'bp_macs2', 'bp_genrich', 'bp_overlap']

        number_peaks.to_csv('evaluate_peaks/Examine_peaks_results_genrich_macs2_union.txt', sep='\t', index = False)



def compare_macs2_combined():
        median = []
        number_peaks = []

        macs_peaks = pd.read_csv('%s/union-peaks.bed' % (macs_dir), sep='\t', header = None)
        combined_peaks = pd.read_csv('%s/combined/union-peaks_combined.bed' % (macs_dir), sep='\t', header = None)
        median.append([np.median(np.array(macs_peaks[2] - macs_peaks[1])), np.median(np.array(combined_peaks[2] - combined_peaks[1]))])

        # look at overlap
        macs2_in_combined = pd.read_csv('%s/macs2_in_combined.bed' % (overlap_dir), sep='\t', header = None)
        macs2_in_combined = macs2_in_combined.drop_duplicates()
        combined_in_macs2 = pd.read_csv('%s/combined_in_macs2.bed' % (overlap_dir), sep='\t', header = None)
        combined_in_macs2 = combined_in_macs2.drop_duplicates()

        # overlapping basepairs
        macs2_basepair = pd.read_csv('%s/macs2_peak_regions.bed' % (overlap_dir), sep='\t', header = None)
        combined_basepair = pd.read_csv('%s/combined_peak_regions.bed' % (overlap_dir), sep='\t', header = None)
        overlap_bp = pd.read_csv('%s/macs2_combined_overlap.bed' % (overlap_dir), sep='\t', header = None)

        number_peaks.append([len(macs_peaks), len(combined_peaks), len(set(macs2_in_combined[3])), len(set(combined_in_macs2[4])), np.sum(macs2_basepair[2] - macs2_basepair[1]), np.sum(combined_basepair[2] - combined_basepair[1]), np.sum(overlap_bp[2] - overlap_bp[1])])
        median = np.array(median)
	print(median)

        number_peaks = pd.DataFrame(number_peaks)
        number_peaks.columns = ['N_peaks_macs2', 'N_peaks_combined', 'N_peaks_macs2_in_combined', 'N_peaks_combined_in_macs2', 'bp_macs2', 'bp_combined', 'bp_overlap']

        number_peaks.to_csv('evaluate_peaks/Examine_peaks_results_combined_macs2.txt', sep='\t', index = False)




def compare_genrich_combined():
        median = []
        number_peaks = []

        genrich_peaks = pd.read_csv('%s/union-peaks.bed' % (genrich_dir), sep='\t', header = None)
        combined_peaks = pd.read_csv('%s/combined/union-peaks_combined.bed' % (genrich_dir), sep='\t', header = None)
        median.append([np.median(np.array(genrich_peaks[2] - genrich_peaks[1])), np.median(np.array(combined_peaks[2] - combined_peaks[1]))])

        # look at overlap
        genrich_in_combined = pd.read_csv('%s/genrich_in_combined.bed' % (overlap_dir), sep='\t', header = None)
        genrich_in_combined = genrich_in_combined.drop_duplicates()
        combined_in_genrich = pd.read_csv('%s/combined_in_genrich.bed' % (overlap_dir), sep='\t', header = None)
        combined_in_genrich = combined_in_genrich.drop_duplicates()

        # overlapping basepairs
        genrich_basepair = pd.read_csv('%s/genrich_peak_regions.bed' % (overlap_dir), sep='\t', header = None)
        combined_basepair = pd.read_csv('%s/combined_peak_regions.bed' % (overlap_dir), sep='\t', header = None)
        overlap_bp = pd.read_csv('%s/genrich_combined_overlap.bed' % (overlap_dir), sep='\t', header = None)

        number_peaks.append([len(genrich_peaks), len(combined_peaks), len(set(genrich_in_combined[3])), len(set(combined_in_genrich[3])), np.sum(genrich_basepair[2] - genrich_basepair[1]), np.sum(combined_basepair[2] - combined_basepair[1]), np.sum(overlap_bp[2] - overlap_bp[1])])
        median = np.array(median)
        print(median)

        number_peaks = pd.DataFrame(number_peaks)
        number_peaks.columns = ['N_peaks_genrich', 'N_peaks_combined', 'N_peaks_genrich_in_combined', 'N_peaks_combined_in_genrich', 'bp_genrich', 'bp_combined', 'bp_overlap']

        number_peaks.to_csv('evaluate_peaks/Examine_peaks_results_combined_genrich.txt', sep='\t', index = False)





if __name__ == '__main__':
        macs_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks'
        genrich_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_Genrich'
        overlap_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_Examine'

	#compare_macs2_combined()
	#compare_macs2_genrich_union()
	compare_genrich_combined()





