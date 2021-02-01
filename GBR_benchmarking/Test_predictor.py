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


def test_predictor(sample):

    try:
        test_dat = pd.read_csv('%s/input/%s_training_dat.txt' % (training_model_dir, sample), sep = '\t', index_col = 0)
    except:
        return

    X = test_dat[['Dosage', 'Imputed_dosage', 'DP', 'GQ', 'Distance_to_peaks']]
    X = X.copy()
    X['DP'] = np.power(10, X['DP']) - 1
    X['Inter1'] = X['Dosage'] * X['DP']
    X['Inter2'] = X['Dosage'] * X['GQ']
    X['Inter3'] = X['Imputed_dosage'] * X['DP']
    X['Inter4'] = X['Imputed_dosage'] * X['GQ']


    for model in ['random_forest', 'linear_regression', 'logistic_regression', 'ordinal_regression']:
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

    test_dat.to_csv('%s/output/%s_minDP%d_variants_THR%.1f.txt' % (training_model_dir, sample, minDP, THR), sep='\t')
    



if __name__ == '__main__':

    minDP = int(sys.argv[1])
    THR = float(sys.argv[2])
    training_model_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model/minDP%d' % minDP
    VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/minDP%d' % minDP
    Impute_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation/minDP%d' % minDP
    suffix = '.filtered.minDP%d.imputed.dosage_genotype.bed' % minDP
    save_dir = 'Performance'
    os.makedirs(training_model_dir, exist_ok = True)
    os.makedirs(os.path.join(training_model_dir, 'output'), exist_ok = True)

    test_samples = pd.read_csv('Test_samples.txt', sep='\t', header =None)
    test_samples = np.array(test_samples[0])

    for s in test_samples:
        try:
            dat = pd.read_csv('%s/output/%s_minDP%d_variants_THR%.1f.txt' % (training_model_dir, s, minDP, THR), sep='\t', nrows=10)
        except:
            print('Test for %s' % s)
            test_predictor(s)
