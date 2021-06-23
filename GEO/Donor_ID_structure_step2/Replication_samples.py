import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import pdb
from collections import OrderedDict
import sys

minDP = int(sys.argv[1])
study = sys.argv[2]

GEO_header = pd.read_csv('/work-zfs/abattle4/heyuan/Variant_calling/GEO/SraRunTables/batch1/SraRunTable_Haematopoiesis.csv', sep='\t', index_col = 0)
GEO_info = pd.read_csv('/work-zfs/abattle4/heyuan/Variant_calling/GEO/SraRunTables/SraRunTable_noSC.csv', sep=',', index_col = 0, header = None)
GEO_info.columns = list(GEO_header.columns) + ['NA']
GEO_info = GEO_info.fillna('None')

GEO_info['samples'] = GEO_info.index
GEO_info = GEO_info.set_index('BioProject').loc[study]

#GEO_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/merged_vcf_files'
GEO_dat = pd.DataFrame()
VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/ATAC_seq/alignment_bowtie/VCF_files/minDP%d' % minDP
for chromosome in range(1, 23):
    print(chromosome)
    #GEO_file = '%s/chr%d.maf005.biallelic.recode.dosage.txt' % (GEO_dir, chromosome)
    GEO_file = '%s/dosage_by_sample_matrix_chr%s.recode.txt' % (VCF_dir, chromosome)
    token = pd.read_csv(GEO_file, sep=' ', index_col = 0)
    GEO_dat = GEO_dat.append(token)

samples = [x for x in GEO_dat.columns if x.startswith('SRR')]
samples = np.intersect1d(samples, GEO_info['samples'])
GEO_dat = GEO_dat[samples]

#original_samples = np.array([x[:10] for x in samples])
#GEO_info = GEO_info.loc[original_samples]

GEO_info = GEO_info.set_index('samples').loc[samples]

def compute_spearman_cor_df(dat, row_names):
    dat[dat == -1] = np.nan

    #spearman_cor = spearmanr(dat.values, axis=0)[0]
    #spearman_cor = pd.DataFrame(spearman_cor, columns=row_names, index = row_names)
    #spearman_cor['sample1'] = spearman_cor.index
    #spearman_cor = spearman_cor.melt(id_vars = 'sample1')
    #spearman_cor.columns = ['sample1', 'sample2', 'Spearman_correlation']

    pearson_cor = dat.corr()
    pearson_cor['sample1'] = pearson_cor.index
    pearson_cor = pearson_cor.melt(id_vars = 'sample1')
    pearson_cor.columns = ['sample1', 'sample2', 'Pearson_correlation']

    spearman_cor = pearson_cor.copy()
    spearman_cor.columns = ['sample1', 'sample2', 'Spearman_correlation']

    cor_dfi = spearman_cor.merge(pearson_cor, left_on = ['sample1', 'sample2'], right_on = ['sample1', 'sample2'])
    cor_dfi = cor_dfi[cor_dfi['sample1'] != cor_dfi['sample2']]

    n_variants = pd.DataFrame(dat.apply(lambda x: len(x.dropna())))
    n_variants.columns = ['n_variants']
    cor_dfi = cor_dfi.merge(n_variants, left_on = ['sample1'], right_index = True)

    return cor_dfi

cor_df = compute_spearman_cor_df(GEO_dat, GEO_dat.columns)

spearman_cor = cor_df.pivot(index = 'sample1', columns = 'sample2')['Spearman_correlation']
pearson_cor = cor_df.pivot(index = 'sample1', columns = 'sample2')['Pearson_correlation']
cor_df.to_csv('correlation/correlation_minDP%d_octopus_%s.csv' % (minDP, study), sep='\t')
spearman_cor.to_csv('correlation/correlation_Spearman_minDP%d_octopus_%s.csv' % (minDP, study), sep='\t')
pearson_cor.to_csv('correlation/correlation_Pearson_minDP%d_octopus_%s.csv' % (minDP, study),  sep='\t')

