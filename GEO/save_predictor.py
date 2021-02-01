import os
import sys
import numpy as np
import pandas as pd
from collections import Counter
from sklearn import tree
from sklearn import ensemble
from sklearn import linear_model
import statsmodels.api as sm
import pickle

from utils import *

sys.path.append('/home-net/home-4/yhe23@jhu.edu/miniconda3/lib/python3.6/site-packages')
import graphviz 


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




def save_train_dataset(sample, minDP):

    print('Save training data for %s' % sample)

    golden_standard_genotype = readin_golden_standard_genotype(golden_standard_dir, sample)
    if golden_standard_genotype is None:
        return

    dat_sample = readin_variants(minDP, sample)
    dat_sample = dat_sample.merge(golden_standard_genotype, left_index = True, right_index = True)
    dat_sample['Distance_to_peaks'] = np.log10(dat_sample['Distance_to_peaks'] + 1)
    dat_sample = dat_sample.drop_duplicates()

    dat_sample.to_csv('%s/input/%s_training_dat.txt' % (training_data_dir, sample), sep='\t')
    print('done\n')

    return dat_sample



def test_predictor(sample):

    try:
        test_dat = pd.read_csv('%s/input/%s_training_dat.txt' % (training_data_dir, sample), sep = '\t', index_col = 0)
    except:
        return

    X = test_dat[['Dosage', 'Imputed_dosage', 'DP', 'GQ', 'Distance_to_peaks']]
    X = X.copy()
    X['DP'] = np.power(10, X['DP']) - 1
    X['Inter1'] = X['Dosage'] * X['DP']
    X['Inter2'] = X['Dosage'] * X['GQ']
    X['Inter3'] = X['Imputed_dosage'] * X['DP']
    X['Inter4'] = X['Imputed_dosage'] * X['GQ']


    for model in ['random_forest']:
        clf = pickle.load(open('%s/%s_THR%.1f.model' % (training_model_dir, model, THR) , 'rb'))
        if (model == 'random_forest') | (model == 'logistic_regression') | (model == 'random_forest_maxDP10'):
            Y_prob = clf.predict_proba(X)
            Y_dosage = Y_prob[:, 1]+ Y_prob[:, 2] * 2
            test_dat['Y_%s' % model] = Y_dosage
        elif (model == 'ordinal_regression'):
            test_dat['Y_%s' % model] = clf.predict(X)
        elif (model == 'linear_regression') :
            #X2 = sm.add_constant(X)
            test_dat['Y_%s' % model] = np.array(clf.predict(X))

    test_dat.to_csv('%s/output/%s_minDP%d_variants_THR%.1f.txt' % (training_data_dir, sample, minDP, THR), sep='\t')
    



if __name__ == '__main__':

    minDP = int(sys.argv[1])
    THR = float(sys.argv[2])
    training_model_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model/minDP%d' % minDP

    training_data_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/PRJNA603385/ATAC_seq/alignment_bowtie'

    samples = pd.read_csv('samples.txt', sep='\t', header =None)
    samples = np.array(samples[0])

    for sample in samples:

        try:
            dat_f = pd.read_csv('%s/input/%s_training_dat.txt' % (training_data_dir, sample), sep = '\t', index_col = 0, nrows = 10)
        except:
            print('Save dataframe for %s' % sample)
            dat_f = save_train_dataset(sample, minDP)
            continue

        try:
            dat = pd.read_csv('%s/output/%s_minDP%d_variants_THR%.1f.txt' % (training_data_dir, sample, minDP, THR), sep='\t', nrows=10)
        except:
            print('Predict for %s' % sample)
            test_predictor(sample)
