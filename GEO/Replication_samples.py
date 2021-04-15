import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import pdb

GEO_info = pd.read_csv('/work-zfs/abattle4/heyuan/Variant_calling/GEO/SraRunTables/SraRunTable_Haematopoiesis.csv', 
                       sep='\t', index_col = 0)
GEO_info['Donor'] = ['Donor: %s' % x for x in GEO_info['Donor']]
GEO_info = GEO_info.fillna('None') 


GEO_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/ATAC_seq/alignment_bowtie/Imputation/minDP3'

GBR_dat = pd.DataFrame()
GEO_dat = pd.DataFrame()
for chromosome in range(1, 23):
    print(chromosome)
    GEO_file = '%s/dosage_by_sample_matrix_chr%d.txt' % (GEO_dir, chromosome)
    token = pd.read_csv(GEO_file, sep=' ')
    token = token[token.columns[3:-1]]
    GEO_dat = GEO_dat.append(token)


samples = np.intersect1d(GEO_info.index, GEO_dat.columns)
GEO_dat = GEO_dat[samples]
GEO_info = GEO_info.loc[samples]
nametag = GEO_info[['Cell_type']].apply(lambda x: '; '.join(x), axis=1)

def compute_spearman_cor_df(dat, row_names):
    spearman_cor = spearmanr(dat.values, axis=0)[0]
    spearman_cor = pd.DataFrame(spearman_cor, columns=row_names, index = row_names)
    cor_array = np.array(spearman_cor)[np.tril_indices_from(spearman_cor, k=-1)]

    pearson_cor = dat.corr()
    pearson_cor.columns=row_names
    pearson_cor.index = row_names

    return [spearman_cor, pearson_cor, cor_array]


[spearman_cor, pearson_cor, GEO_cor] = compute_spearman_cor_df(GEO_dat, samples)
spearman_cor.to_csv('SraRunTables/correlation.csv', sep='\t')
pearson_cor.to_csv('SraRunTables/correlation_pearson.csv', sep='\t')

