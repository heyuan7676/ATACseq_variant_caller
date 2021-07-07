import pandas as pd
import numpy as np
import os
import sys
import pdb


if __name__ == '__main__':
    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/ATAC_seq/alignment_bowtie'
    minDP = int(sys.argv[1])
    minMQ = int(sys.argv[2])
    sample = sys.argv[3]

    VCF_dir = os.path.join(root_dir, 'VCF_files', 'minDP%d' % minDP)
    output_file = '%s/%s.filtered.minDP%d.recode.dosage_genotype.octopus.bed' % (VCF_dir, sample, minDP)
    if os.path.isfile(output_file):
        sys.exit("Exist")

    try:
        pp_dat = pd.read_csv('%s/VCF_files/%s.octopus.fast.vcf' % (root_dir, sample), sep='\t', header = None,  comment = "#")
    except:
        print('%s/VCF_files/%s.octopus.fast.vcf' % (root_dir, sample))
        sys.exit("Octopus not available")

    gt_vector = np.array([(x.split(':')[0]) for x in np.array(pp_dat[9])])
    valid_loci = np.where(np.array(['.' not in x for x in gt_vector]))[0]
    pp_dat = pp_dat.iloc[valid_loci]
    gt_vector = np.array([(x.split(':')[0]) for x in np.array(pp_dat[9])])
    pp_dat['GT'] = [int(x.split('|')[0]) + int(x.split('|')[1]) for x in gt_vector]
    pp_dat['GQ'] = [int(x.split(':')[1]) for x in np.array(pp_dat[9])]
    pp_dat['DP'] = [int(x.split(':')[2]) for x in np.array(pp_dat[9])]
    pp_dat['MQ'] = [int(x.split(':')[3]) for x in np.array(pp_dat[9])]

    # filter on minDP
    pp_dat_use = pp_dat[pp_dat['DP'] >= minDP]

    dosage_with_ll = pp_dat_use[['GT', 'GT']]
    dosage_with_ll.columns = ['Dosage', 'GT']
    dosage_with_ll['CHR_POS'] = ['_'.join(map(str, x)).replace('chr', '') for x in zip(np.array(pp_dat_use[0]), np.array(pp_dat_use[1]))]
    dosage_with_ll = dosage_with_ll[['CHR_POS', 'Dosage', 'GT']]

    VCF_dir = os.path.join(root_dir, 'VCF_files', 'minDP%d' % minDP)
    os.makedirs(VCF_dir, exist_ok = True)
    dosage_with_ll.to_csv('%s/%s.filtered.minDP%d.recode.dosage_genotype.octopus.bed' % (VCF_dir, sample, minDP), sep='\t', index = False)


    # filter on MQ
    pp_dat = pp_dat[pp_dat['MQ'] >= minMQ]

    dosage_with_ll = pp_dat[['GT', 'GT']]
    dosage_with_ll.columns = ['Dosage', 'GT']
    dosage_with_ll['CHR_POS'] = ['_'.join(map(str, x)).replace('chr', '') for x in zip(np.array(pp_dat[0]), np.array(pp_dat[1]))]
    dosage_with_ll = dosage_with_ll[['CHR_POS', 'Dosage', 'GT']]

    VCF_dir = os.path.join(root_dir, 'VCF_files', 'minMQ%d' % minMQ)
    os.makedirs(VCF_dir, exist_ok = True)
    dosage_with_ll.to_csv('%s/%s.filtered.minMQ%d.recode.dosage_genotype.octopus.bed' % (VCF_dir, sample, minMQ), sep='\t', index = False)



