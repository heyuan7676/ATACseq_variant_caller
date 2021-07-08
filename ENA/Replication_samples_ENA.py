import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import pdb
from collections import OrderedDict
import sys

#minDP = int(sys.argv[1])
#study = sys.argv[2]

minDP = 5
study = 'PRJEB17990'

ENA_info = pd.read_csv('/work-zfs/abattle4/heyuan/Variant_calling/ENA/records/sample_metadata_noSC.tsv', sep='\t', index_col = 0)
ENA_info = ENA_info.fillna('None')

ENA_info['samples'] = ENA_info['run_accession']
ENA_info = ENA_info.loc[study]

VCF_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/ENA/ATAC_seq/alignment_bowtie/VCF_files/minDP%d' % minDP
ENA_dat = pd.DataFrame()
samples_vec = []
for sample in ENA_info['samples']:
    ENA_file = '%s/%s.filtered.minDP%d.recode.dosage_genotype.octopus.bed' % (VCF_dir, sample, minDP)
    try:
        token = pd.read_csv(ENA_file, sep='\t', index_col = 0)
        token = token[['Dosage']]
    except:
        print('%s has no vcf file' % sample)
        continue
    if len(ENA_dat) == 0:
        ENA_dat = token.copy()
    else:
        ENA_dat = ENA_dat.merge(token, left_index = True, right_index = True)
    samples_vec.append(sample)

ENA_dat.columns = samples_vec
ENA_info = ENA_info.set_index('samples').loc[samples_vec]

def compute_spearman_cor_df(dat, row_names):
    dat[dat == -1] = np.nan

    spearman_cor = spearmanr(dat.values, axis=0)[0]
    spearman_cor = pd.DataFrame(spearman_cor, columns=row_names, index = row_names)
    spearman_cor['sample1'] = spearman_cor.index
    spearman_cor = spearman_cor.melt(id_vars = 'sample1')
    spearman_cor.columns = ['sample1', 'sample2', 'Spearman_correlation']

    pearson_cor = dat.corr()
    pearson_cor['sample1'] = pearson_cor.index
    pearson_cor = pearson_cor.melt(id_vars = 'sample1')
    pearson_cor.columns = ['sample1', 'sample2', 'Pearson_correlation']

    cor_dfi = spearman_cor.merge(pearson_cor, left_on = ['sample1', 'sample2'], right_on = ['sample1', 'sample2'])
    cor_dfi = cor_dfi[cor_dfi['sample1'] != cor_dfi['sample2']]

    n_variants = pd.DataFrame(dat.apply(lambda x: len(x.dropna())))
    n_variants.columns = ['n_variants']
    cor_dfi = cor_dfi.merge(n_variants, left_on = ['sample1'], right_index = True)

    return cor_dfi

cor_df = compute_spearman_cor_df(ENA_dat, ENA_dat.columns)

spearman_cor = cor_df.pivot(index = 'sample1', columns = 'sample2')['Spearman_correlation']
pearson_cor = cor_df.pivot(index = 'sample1', columns = 'sample2')['Pearson_correlation']
cor_df.to_csv('correlation/correlation_minDP%d_octopus_%s.csv' % (minDP, study), sep='\t')
spearman_cor.to_csv('correlation/correlation_Spearman_minDP%d_octopus_%s.csv' % (minDP, study), sep='\t')
pearson_cor.to_csv('correlation/correlation_Pearson_minDP%d_octopus_%s.csv' % (minDP, study),  sep='\t')

