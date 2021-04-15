import pandas as pd
import numpy as np
import sys
import pdb

chromosome = int(sys.argv[1])

VCF_filename = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/ATAC_seq/alignment_bowtie/VCF_files/merged_vcf_chromosome%d.vcf' % chromosome

header = []
with open(VCF_filename,'rt') as f:
    for l in f:
        if l.rstrip().startswith('#'):
            header.append(l.rstrip())
        else:
            break

VCF_dat = pd.read_csv(VCF_filename, usecols=[0,1,2,3,4,5,6,7,8], comment='#', header=None, sep='\t')
VCF_dat[2] = VCF_dat[[0,1]].apply(lambda x: '_'.join((str(x[0]), str(x[1]))), axis=1)
VCF_dat[0] = ['chr%d' % x for x in VCF_dat[0]]

# Imputation

dosage_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/ATAC_seq/alignment_bowtie/Imputation/minDP3'
dosage_file = '%s/dosage_by_sample_matrix_chr%d.txt' % (dosage_dir, chromosome)
dosage_dat = pd.read_csv(dosage_file, sep=' ', index_col = [0])
samples = [x for x in dosage_dat.columns if x.startswith('SRR')]

merged_dat = VCF_dat.merge(dosage_dat[samples], left_on=2, right_index=True)
merged_dat[8] = ['DS'] * len(merged_dat)
merged_dat = merged_dat[merged_dat[samples].apply(lambda x: np.sum(x) > 3, axis=1)]


# modify the last line of header
header[-1] = '\t'.join([x for x in header[-1].split('\t') if not x.startswith('SRR')] + list(samples))

output_file = '%s/dosage_by_sample_matrix_chr%d.recode.forFastQTL.vcf' % (dosage_dir, chromosome)
print(output_file)

f = open(output_file, "w")
for l in header:
    f.write(l + '\n')


for l in np.array(merged_dat):
    f.write('\t'.join(map(str, l)) + '\n')

f.close()



