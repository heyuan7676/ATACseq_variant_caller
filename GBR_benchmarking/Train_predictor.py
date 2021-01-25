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




def train_predictor(minDP):
    
    files = os.listdir(Impute_dir)
    files = np.sort([f for f in files if suffix in f])
    
    # devide into training and testing
    np.random.seed(1)
    training_files = np.random.choice(files, int(len(files) * 0.8), replace = False)
    test_files = np.sort(list(set(files) - set(training_files)))
    
    # read in data
    training_dat = pd.DataFrame()
    for f in np.sort(training_files):
        print('Read in data for %s' % f)
        sample = f.replace(suffix, '')
        golden_standard_genotype = readin_golden_standard_genotype(golden_standard_dir, sample)
        if golden_standard_genotype is not None:
            dat_sample = readin_variants(minDP, sample)
            dat_sample = dat_sample[np.abs(dat_sample['Imputed_dosage'] - dat_sample['Dosage']) > THR]
            dat_sample = dat_sample.merge(golden_standard_genotype, left_index = True, right_index = True)
            training_dat = training_dat.append(dat_sample)
        print("\n")

    training_dat['Distance_to_peaks'] = np.log10(training_dat['Distance_to_peaks'] + 1)
    training_dat = training_dat.drop_duplicates()
    
    # train the model
    X = training_dat[['Dosage', 'Imputed_dosage', 'DP', 'GQ', 'Distance_to_peaks']]
    Y = np.array(training_dat['True_GT'])

    # Regression
    print('Fit in linear regression')
    lr = linear_model.LinearRegression().fit(X, Y)
    pickle.dump(lr, open('%s/linear_regression_THR%.1f.model' % (training_model_dir, THR), 'wb'))
    Yhat = lr.predict(X)
    mse_lr = np.sum([x**2 for x in np.array(Yhat - Y)]) / len(Y)

    # logistic Regression
    print('Fit in logistic regression')
    lr = linear_model.LogisticRegression(multi_class = 'multinomial', solver = 'lbfgs').fit(X, Y)
    pickle.dump(lr, open('%s/logistic_regression_THR%.1f.model' % (training_model_dir, THR), 'wb'))
    Yhat = lr.predict(X)
    mse_lgr = np.sum([x**2 for x in np.array(Yhat - Y)]) / len(Y)


    # Random forest
    print('Fit in random forest')
    clf = ensemble.RandomForestClassifier(max_depth = 3).fit(X, Y)
    pickle.dump(clf, open('%s/random_forest_THR%.1f.model' % (training_model_dir, THR), 'wb'))
    Yhat = clf.predict(X)
    mse_clf = np.sum([x**2 for x in np.array(Yhat - Y)]) / len(Y)
    
    print('On training data, MSE = %.2f using linear regression, = %.2f using logistic regression, and = %.2f using random forest' % (mse_lr, mse_lgr, mse_clf))
    #tree.export_graphviz(clf, out_file="%s/tree.dot" % training_model_dir, feature_names = X.columns,  class_names=['0', '1', '2'], filled = True)

    test_samples = [x.replace(suffix, '') for x in test_files]
    for s in test_samples:
        try:
            dat = pd.read_csv('%s/output/%s_minDP%d_variants.txt' % (training_model_dir, s, minDP), sep='\t')
        except:
            print('Test for %s' % s)
	    # golden standard
            golden_standard_genotype = readin_golden_standard_genotype(golden_standard_dir, s)
            golden_variants = pd.DataFrame.from_dict(Counter(golden_standard_genotype['True_GT']), orient='index')

            test_predictor(s, golden_standard_genotype, golden_variants, model = 'random_forest')
            test_predictor(s, golden_standard_genotype, golden_variants, model = 'linear_regression')
            test_predictor(s, golden_standard_genotype, golden_variants, model = 'logistic_regression')


    pd.DataFrame(test_samples).to_csv('Test_samples.txt', index = False, header = False)

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

    test_samples = train_predictor(minDP)
    for s in test_samples:
        try:
            dat = pd.read_csv('%s/output/%s_minDP%d_variants.txt' % (training_model_dir, s, minDP), sep='\t')
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


