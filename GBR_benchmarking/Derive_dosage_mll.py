import pandas as pd
import numpy as np
import os
import sys
import pdb


def achieve_dosage_and_ll(rowi):
    '''
    Map PL scores back to likelihood
    Output: expected genotype / dosage of the variant, and the highest likelihood
    '''
    pl_scores = list(map(int, str(rowi).split(',')))
    pl_scores = np.array(pl_scores)
    if(len(pl_scores)) > 3:
        return [-1, -1]
    
    # avoid overflow
    pl_scores[pl_scores > 20] = 20
    post_pp = np.array(list(map(lambda x: 1.0 / (10 ** x), pl_scores)))
    post_pp = post_pp / float(np.sum(post_pp))

    dosage = post_pp[0] * 0 +post_pp[1] * 1 + post_pp[2] * 2
    if len(set(post_pp)) < len(post_pp):
        pool = np.where(post_pp == np.max(post_pp))[0]
        max_post_pp = np.random.choice(pool, size = 1)[0]
    else:
        max_post_pp = np.argmax(post_pp)

    return [dosage, max_post_pp]


if __name__ == '__main__':
    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie'
    minDP = int(sys.argv[1])
    sample = sys.argv[2]

    VCF_dir = os.path.join(root_dir, 'VCF_files', 'minDP%d' % minDP)
    filename = '%s/%s.filtered.minDP2.recode.vcf.gz' % (VCF_dir, sample)
    pp_dat = pd.read_csv(filename, sep='\t', header = None,  compression='gzip', comment = "#")
    pp_dat['PL_index'] = [x.split(':').index('PL') for x in np.array(pp_dat[8])]
    pp_dat_list = pp_dat.apply(lambda x: x[9].split(':')[x['PL_index']], axis = 1) 

    dosage_with_ll = pd.DataFrame(list(map(achieve_dosage_and_ll, np.array(pp_dat_list))))
    dosage_with_ll.index = pp_dat[0]
    dosage_with_ll.to_csv('%s/%s.filtered.recode.INFO.formatted.variants.dosage.txt' % (VCF_dir, sample), sep='\t')


