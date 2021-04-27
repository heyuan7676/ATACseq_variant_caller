import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import pdb
from collections import OrderedDict

GEO_header = pd.read_csv('/work-zfs/abattle4/heyuan/Variant_calling/GEO/SraRunTables/batch1/SraRunTable_Haematopoiesis.csv', sep='\t', index_col = 0)
GEO_info = pd.read_csv('/work-zfs/abattle4/heyuan/Variant_calling/GEO/SraRunTables/SraRunTable_noSC.csv', sep=',', index_col = 0, header = None)
GEO_info.columns = list(GEO_header.columns) + ['NA']
GEO_info = GEO_info.fillna('None')


GEO_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/merged_vcf_files'
GEO_dat = pd.DataFrame()
header = pd.read_csv('%s/samples.txt' % GEO_dir, header = None)
header = np.array(header[0])
for chromosome in range(22, 23):
    print(chromosome)
    GEO_file = '%s/chr%d.maf005.biallelic.recode.dosage.txt' % (GEO_dir, chromosome)
    token = pd.read_csv(GEO_file, sep='\t', header = None, index_col = 0)
    GEO_dat = GEO_dat.append(token)

GEO_dat.columns = header
samples = np.intersect1d(GEO_info.index, GEO_dat.columns)
GEO_dat = GEO_dat[samples]
GEO_info = GEO_info.loc[samples]

def compute_spearman_cor_df(dat, row_names, study_name):
    spearman_cor = spearmanr(dat.values, axis=0)[0]
    spearman_cor = pd.DataFrame(spearman_cor, columns=row_names, index = row_names)
    spearman_cor['sample1'] = spearman_cor.index
    spearman_cor = spearman_cor.melt(id_vars = 'sample1')
    spearman_cor.columns = ['sample1', 'sample2', 'Spearman_correlation']

    pearson_cor = dat.corr()
    pearson_cor.columns=row_names
    pearson_cor.index = row_names
    pearson_cor['sample1'] = pearson_cor.index
    pearson_cor = pearson_cor.melt(id_vars = 'sample1')
    pearson_cor.columns = ['sample1', 'sample2', 'Pearson_correlation']

    cor_dfi = spearman_cor.merge(pearson_cor, left_on = ['sample1', 'sample2'], right_on = ['sample1', 'sample2'])
    cor_dfi = cor_dfi[cor_dfi['sample1'] != cor_dfi['sample2']]
    cor_dfi['study'] = study_name

    return cor_dfi



cor_df = pd.DataFrame()
for study in GEO_info.groupby('BioProject'):
    cor_dfi = compute_spearman_cor_df(GEO_dat[study[1].index], study[1].index, study[0])
    cor_df = cor_df.append(cor_dfi)

spearman_cor = cor_df.pivot(index = 'sample1', columns = 'sample2')['Spearman_correlation']
pearson_cor = cor_df.pivot(index = 'sample1', columns = 'sample2')['Pearson_correlation']
cor_df.to_csv('SraRunTables/correlation.csv', sep='\t')
spearman_cor.to_csv('SraRunTables/correlation_Spearman.csv', sep='\t')
pearson_cor.to_csv('SraRunTables/correlation_Pearson.csv', sep='\t')

