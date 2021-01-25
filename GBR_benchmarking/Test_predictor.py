import os
import sys
import numpy as np
import pandas as pd
from collections import Counter
from sklearn import tree
from sklearn import ensemble
from sklearn import linear_model
import pickle

from utils import *

sys.path.append('/home-net/home-4/yhe23@jhu.edu/miniconda3/lib/python3.6/site-packages')
import graphviz 




def obtain_test_samples(minDP):
    
    files = os.listdir(Impute_dir)
    files = np.sort([f for f in files if suffix in f])
    
    # devide into training and testing
    np.random.seed(1)
    training_files = np.random.choice(files, int(len(files) * 0.8), replace = False)
    test_files = np.sort(list(set(files) - set(training_files)))
    test_samples = [x.replace(suffix, '') for x in test_files]

    return test_samples



def readin_variants(minDP, sample):

    # genotype caller
    [called_genotype, variants_count] = readin_variant_caller_genotype(VCF_dir, minDP, sample)
    called_genotype = append_features(called_genotype, sample, minDP, imputation = False)

    # imputation
    [imputed_genotype, variants_count_imputed] = readin_imputation_genotype(Impute_dir, minDP, sample)

    variants_genotype = imputed_genotype.merge(called_genotype, left_index = True, right_index = True)
    return variants_genotype 



 
def append_features(dati, sample, minDP, imputation = False):
    if not imputation:
        # read in read depth and PL score
        read_depth = pd.read_csv('%s/%s.filtered.recode.INFO.txt' % (VCF_dir, sample), header = None, sep=' ', index_col = 0)
        read_depth.columns = ['DP', 'GQ']
        dati['DP'] = read_depth.loc[np.array(dati.index)]['DP']
        dati['GQ'] = read_depth.loc[np.array(dati.index)]['GQ']

        # read in variants' distance to peaks
        distance_to_peak = pd.read_csv('%s/%s.filtered.minDP%d.recode.variants.toMACS2Peaks.txt' % (VCF_dir, sample, minDP), header = None, sep=' ', index_col = 0, usecols = [0, 2])
        distance_to_peak.columns = ['Distance_to_peaks']
        min_distance = distance_to_peak.groupby(distance_to_peak.index).min()
        min_distance = min_distance.reset_index().drop_duplicates().set_index(0)
        dati['Distance_to_peaks'] = np.array(min_distance.loc[np.array(dati.index)]['Distance_to_peaks'])

    elif imputation:
        distance_to_peak = pd.read_csv('%s/%s.filtered.minDP%d.imputed.variants.toMACS2Peaks.txt' % (Impute_dir, sample, minDP), header = None, sep=' ', index_col = 0, usecols = [0, 2])
        distance_to_peak.columns = ['Distance_to_peaks']
        min_distance_impute = distance_to_peak.groupby(distance_to_peak.index).min()
        min_distance_impute  = min_distance_impute.reset_index().drop_duplicates().set_index(0)
        dati['Distance_to_peaks'] = np.array(min_distance_impute.loc[np.array(dati.index)]['Distance_to_peaks'])


    return dati



 
def test_predictor(sample, golden_standard_genotype, golden_variants, model = 'random_forest'):

    test_dat = readin_variants(minDP, sample)
    test_dat = test_dat.merge(golden_standard_genotype, left_index = True, right_index = True)
    test_dat['Distance_to_peaks'] = np.log10(test_dat['Distance_to_peaks'] + 1)

    X = test_dat[['Dosage', 'Imputed_dosage', 'DP', 'GQ', 'Distance_to_peaks']]

    for model in ['random_forest', 'linear_regression', 'logistic_regression']:
        clf = pickle.load(open('%s/%s_THR%.1f.model' % (training_model_dir, model, THR) , 'rb'))
        if (model == 'random_forest') | (model == 'logistic_regression'):
            Y_prob = clf.predict_proba(X)
            Y_dosage = Y_prob[:, 1]+ Y_prob[:, 2] * 2
            test_dat['Y_%s' % model] = Y_dosage
        elif model == 'linear_regression':
            test_dat['Y_%s' % model] = clf.predict(X)

    test_dat.to_csv('%s/output/%s_minDP%d_variants_THR%.1f.txt' % (training_model_dir, sample, minDP, THR), sep='\t')
    



if __name__ == '__main__':

    minDP = int(sys.argv[1])
    THR = float(sys.argv[2])
    golden_standard_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'
    training_model_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model/minDP%d' % minDP
    VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/minDP%d' % minDP
    Impute_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation/minDP%d' % minDP
    suffix = '.filtered.minDP%d.imputed.dosage_genotype.bed' % minDP
    save_dir = 'Performance'
    os.makedirs(training_model_dir, exist_ok = True)
    os.makedirs(os.path.join(training_model_dir, 'output'), exist_ok = True)

    test_samples = obtain_test_samples(minDP)
    print(test_samples)
    for s in test_samples:
        try:
            dat = pd.read_csv('%s/output/%s_minDP%d_variants_THR%.1f.txt' % (training_model_dir, s, minDP, THR), sep='\t', nrows=10)
        except:
            print('Test for %s' % s)
	    # golden standard
            golden_standard_genotype = readin_golden_standard_genotype(golden_standard_dir, s)
            if golden_standard_genotype is None:
                continue
            golden_variants = pd.DataFrame.from_dict(Counter(golden_standard_genotype['True_GT']), orient='index')

            test_predictor(s, golden_standard_genotype, golden_variants, model = 'random_forest')
            test_predictor(s, golden_standard_genotype, golden_variants, model = 'linear_regression')
            test_predictor(s, golden_standard_genotype, golden_variants, model = 'logistic_regression')


