import pdb
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ranksums
from Evaluation_metrics import *

def retrive_coverage_for_all_samples_imputed(restrict_to_SNP, saveFile):
    samples = pd.read_csv('samples.txt', header=None)
    samples = np.array(samples[0])

    performance = []
    for sample in samples[:1]:
	print('For %s:' % sample)
        # WGS genotype
        WGS_df = read_in_WGS_GT(sample)

        # variant calling information
        #try:
	if 1:
	    print('Evaluate genotype originally called from ATAC-seq reads...')
            orginally_called = obtain_atac_variants_df(sample, WGS_result=WGS_df, restrict_to_SNP=restrict_to_SNP, Imputed = False, return_metric = False)

	    print('Evalutae genotype imputed')
	    imputed = obtain_atac_variants_df(sample, WGS_result=WGS_df, restrict_to_SNP=restrict_to_SNP, Imputed = True, return_metric = False)

	    pdb.set_trace()

	    performance_original = compute_metric(intersection_SNPs, sample = sample, minDP = minDP, Imputed = False)
	    performance_imputed = compute_metric(intersection_SNPs, sample = sample, minDP = 2, Imputed = True)

        #except:
        #    print('%s does not have genotype data' % s)

    performance = pd.DataFrame(performance)
    performance.columns = ['Sample', 'minDP',  
			   'N_variants_by_WGS', 'N_variants_by_ATAC', 'N_overlap_variants_by_ATAC', 'Recovered_percentage', 'Correctly_performance_percentage',
 			   'Recall_AA', 'Recall_AB', 'Recall_BB', 'Precision_AA', 'Precision_AB', 'Precision_BB']
    performance.to_csv(saveFile, sep='\t', index = False)

    return performance




if __name__ == '__main__':
    minDP = 2
    retrive_coverage_for_all_samples_imputed(restrict_to_SNP = True, saveFile = 'performance/All_SNPs_imputed.txt')
