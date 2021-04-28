import numpy as np
import pandas as pd
import sys

sample_structure = pd.read_csv('/work-zfs/abattle4/heyuan/Variant_calling/GEO/SraRunTables/correlation_Spearman_sample_structure.csv', sep='\t')

samples = pd.read_csv('samples.txt', header = None)
samples = np.array(samples[0])

chromosome = int(sys.argv[1])
dosage = pd.read_csv('chr%d.maf005.biallelic.recode.dosage.txt' % chromosome, 
                     header = None, sep='\t', nrows = 10, index_col = 0)
dosage.columns = samples

dosage = dosage.transpose()
dosage = dosage.copy()

dosage['Donor'] = np.array(sample_structure.loc[samples]['sample_levels'])
mean_dosage = dosage.groupby('Donor').mean()
mean_dosage = mean_dosage.transpose()
mean_dosage.to_csv('chr%d.maf005.biallelic.recode.dosage.donors.txt' % chromosome, sep='\t')

