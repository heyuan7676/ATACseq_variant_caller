import os
import sys
sys.path.append('/home-net/home-4/yhe23@jhu.edu/miniconda3/lib/python3.6/site-packages')

import numpy as np
import pandas as pd
from collections import Counter
from sklearn import tree
from sklearn import ensemble
from sklearn import linear_model
import mord as m
import pickle

from utils import *

import graphviz 
import statsmodels.api as sm




def readin_variants(minDP, sample):

    # genotype caller
    [called_genotype, variants_count] = readin_variant_caller_genotype(VCF_dir, minDP, sample)
    called_genotype = append_features(called_genotype, sample, minDP, imputation = False)

    # imputation
    [imputed_genotype, variants_count_imputed] = readin_imputation_genotype(Impute_dir, minDP, sample)

    variants_genotype = imputed_genotype.merge(called_genotype, left_index = True, right_index = True)
    return variants_genotype




def save_train_dataset(sample, minDP):
    
    print('Save training data for %s' % sample)

    golden_standard_genotype = readin_golden_standard_genotype(golden_standard_dir, sample)
    if golden_standard_genotype is None:
        return

    dat_sample = readin_variants(minDP, sample)
    dat_sample = dat_sample.merge(golden_standard_genotype, left_index = True, right_index = True)
    dat_sample['Distance_to_peaks'] = np.log10(dat_sample['Distance_to_peaks'] + 1)
    dat_sample = dat_sample.drop_duplicates()

    dat_sample.to_csv('%s/input/%s_training_dat.txt' % (training_model_dir, sample), sep='\t')
    print('done\n')

    return dat_sample



def train_predictor(minDP):

    samples = pd.read_csv('samples.txt', sep='\t', header= None)
    samples = np.array(samples[0]) 

    # prepare the data
    for sample in np.sort(samples):
        # Read in
        try:
            dat_f = pd.read_csv('%s/input/%s_training_dat.txt' % (training_model_dir, sample), sep = '\t', index_col = 0, nrows = 10)
        except:
            print('Save dataframe for %s' % sample)
            dat_f = save_train_dataset(sample, minDP)
            continue

    # devide into training and testing
    np.random.seed(1)
    training_samples = np.random.choice(samples, int(len(samples) * 0.8), replace = False)
    test_samples = np.sort(list(set(samples) - set(training_samples)))

    training_dat = pd.DataFrame()
    for sample in np.sort(training_samples):
        try:
            dat_f = pd.read_csv('%s/input/%s_training_dat.txt' % (training_model_dir, sample), sep = '\t', index_col = 0)
            dat_f = dat_f[np.abs(dat_f['Dosage'] - dat_f['Imputed_dosage']) >= THR]
            training_dat = training_dat.append(dat_f)
        except:
            print('Training data for %d does not exist' % sample)
            continue

    # train the model
    X = training_dat[['Dosage', 'Imputed_dosage', 'DP', 'GQ', 'Distance_to_peaks']]
    X = X.copy()
    X['DP'] = np.power(10, X['DP']) - 1
    X['Inter1'] = X['Dosage'] * X['DP']
    X['Inter2'] = X['Dosage'] * X['GQ']
    X['Inter3'] = X['Imputed_dosage'] * X['DP']
    X['Inter4'] = X['Imputed_dosage'] * X['GQ']

    Y = np.array(training_dat['True_GT'])

    # Regression
    print('Fit in linear regression')
    lr = linear_model.LinearRegression().fit(X, Y)
    pickle.dump(lr, open('%s/linear_regression_THR%.1f.model' % (training_model_dir, THR), 'wb'))
    Yhat = lr.predict(X)
    #X2 = sm.add_constant(X)
    #est = sm.OLS(Y, X2)
    #est2 = est.fit()
    #pickle.dump(est2, open('%s/linear_regression_THR%.1f.model' % (training_model_dir, THR), 'wb'))
    #Yhat = est2.predict(X2)
    mse_lr = np.sum([x**2 for x in np.array(Yhat - Y)]) / len(Y)


    # logistic Regression
    print('Fit in logistic regression')
    lr = linear_model.LogisticRegression() #Default parameters: alpha=1.0, verbose=0, maxiter=10000
    lr.fit(X, Y)
    pickle.dump(lr, open('%s/logistic_regression_THR%.1f.model' % (training_model_dir, THR), 'wb'))
    Yhat = lr.predict(X)
    mse_lgr = np.sum([x**2 for x in np.array(Yhat - Y)]) / len(Y)


    # logistic Regression
    print('Fit in ordinal regression')
    lr = m.OrdinalRidge() #Default parameters: alpha=1.0, verbose=0, maxiter=10000
    lr.fit(X, Y)
    pickle.dump(lr, open('%s/ordinal_regression_THR%.1f.model' % (training_model_dir, THR), 'wb'))
    Yhat = lr.predict(X)
    mse_or = np.sum([x**2 for x in np.array(Yhat - Y)]) / len(Y)

    # Random forest
    print('Fit in random forest')
    clf = ensemble.RandomForestClassifier(max_depth = 10).fit(X, Y)
    pickle.dump(clf, open('%s/random_forest_THR%.1f.model' % (training_model_dir, THR), 'wb'))
    Yhat = clf.predict(X)
    mse_clf = np.sum([x**2 for x in np.array(Yhat - Y)]) / len(Y)
    
    print('On training data, MSE = %.2f using linear regression, = %.2f using logistic regression, and = %.2f using ordinal regression, and = %.2f using random forest' % (mse_lr, mse_lgr, mse_or, mse_clf))
    pd.DataFrame(test_samples).to_csv('Test_samples.txt', index = False, header = False)




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



 
if __name__ == '__main__':

    minDP = int(sys.argv[1])
    THR = float(sys.argv[2])
    golden_standard_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'

    training_model_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Model/minDP%d' % minDP
    VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/minDP%d' % minDP
    Impute_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation/minDP%d' % minDP

    save_dir = 'Performance'
    os.makedirs(training_model_dir, exist_ok = True)
    os.makedirs(os.path.join(training_model_dir, 'input'), exist_ok = True)
    os.makedirs(os.path.join(training_model_dir, 'output'), exist_ok = True)

    train_predictor(minDP)
