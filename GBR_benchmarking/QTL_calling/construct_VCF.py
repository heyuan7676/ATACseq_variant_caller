import pandas as pd
import numpy as np
import gzip
import sys
import pdb

chromosome = int(sys.argv[1])

VCF_filename = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.vcf.withchr.gz' % chromosome

header = []
with gzip.open(VCF_filename,'rt') as f:
    for l in f:
        if l.rstrip().startswith('#'):
            header.append(l.rstrip())
        else:
            break
samples = [x for x in header[-1].split('\t') if x.startswith('HG')]

VCF_dat = pd.read_csv(VCF_filename, usecols=[0,1,2,3,4,5,6,7,8], comment='#', header=None, sep='\t')
VCF_dat.head()


# Genotype Caller

dosage_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/minDP3'
dosage_file = '%s/dosage_by_sample_matrix_chr%d.recode.txt' % (dosage_dir, chromosome)
dosage_dat = pd.read_csv(dosage_file, sep=' ', index_col = [0])


merged_dat = VCF_dat.merge(dosage_dat[samples], left_on=2, right_index=True)
merged_dat[8] = ['DS'] * len(merged_dat)
merged_dat = merged_dat.replace(-1, '.')
merged_dat = merged_dat[merged_dat[samples].apply(lambda x: np.sum(x!='.') > 50, axis=1)]
merged_dat = merged_dat[merged_dat[samples].apply(lambda x: np.sum(x[x!='.']) > 3, axis=1)]


output_file = '%s/dosage_by_sample_matrix_chr%d.recode.forFastQTL.vcf' % (dosage_dir, chromosome)
print(output_file)

f = open(output_file, "w")
for l in header:
    f.write(l + '\n')
    
    
for l in np.array(merged_dat):
    f.write('\t'.join(map(str, l)) + '\n')
    
f.close()


# Imputation

dosage_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation/minDP3'
dosage_file = '%s/dosage_by_sample_matrix_chr%d.txt' % (dosage_dir, chromosome)
dosage_dat = pd.read_csv(dosage_file, sep=' ', index_col = [0])


merged_dat = VCF_dat.merge(dosage_dat[samples], left_on=2, right_index=True)
merged_dat[8] = ['DS'] * len(merged_dat)
merged_dat = merged_dat[merged_dat[samples].apply(lambda x: np.sum(x) > 3, axis=1)]


output_file = '%s/dosage_by_sample_matrix_chr%d.recode.forFastQTL.vcf' % (dosage_dir, chromosome)
print(output_file)

f = open(output_file, "w")
for l in header:
    f.write(l + '\n')


for l in np.array(merged_dat):
    f.write('\t'.join(map(str, l)) + '\n')

f.close()



# Integration
dosage_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Integration/minDP3'
dosage_file = '%s/dosage_by_sample_matrix_chr%d.THR0.0_Y_random_forest.txt' % (dosage_dir, chromosome)
dosage_dat = pd.read_csv(dosage_file, sep=' ', index_col = [0])


merged_dat = VCF_dat.merge(dosage_dat[samples], left_on=2, right_index=True)
merged_dat[8] = ['DS'] * len(merged_dat)
merged_dat = merged_dat.replace(-1, '.')
merged_dat = merged_dat[merged_dat[samples].apply(lambda x: np.sum(x!='.') > 50, axis=1)]
merged_dat = merged_dat[merged_dat[samples].apply(lambda x: np.sum(x[x!='.']) > 3, axis=1)]


output_file = '%s/dosage_by_sample_matrix_chr%d.recode.forFastQTL.vcf' % (dosage_dir, chromosome)
print(output_file)

f = open(output_file, "w")
for l in header:
    f.write(l + '\n')


for l in np.array(merged_dat):
    f.write('\t'.join(map(str, l)) + '\n')

f.close()

