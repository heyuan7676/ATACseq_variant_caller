import sys
import numpy as np
import pandas as pd

sys.path.append('/work-zfs/abattle4/heyuan/Variant_calling/scripts/QTL')
from QTL_calling import compute_QTL_peakLevel

if __name__ == '__main__':
    ROOT_DIR = sys.argv[1]  # /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/
    peak_calling = sys.argv[2]
    WINDOW = int(sys.argv[3])
    minDP = int(sys.argv[4])

    alignment_dir = 'alignment_bowtie' # alignment_subsample_0.5
    GT_subDir = 'minDP%d' % minDP

    root_dir = '%s/%s' % (ROOT_DIR, alignment_dir)
    VCF_dir = '%s/VCF_files' % root_dir
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'

    if peak_calling == 'macs2':
        PEAK_dir = '%s/Peaks_MACS2' % root_dir
        QTL_dir = '%s/QTLs/%s' % (root_dir, GT_subDir)

    if peak_calling == 'macs2_combined':
        PEAK_dir = '%s/Peak_version1/Peaks/combined' % root_dir
        QTL_dir = '%s/QTLs_combined/%s' % (root_dir, GT_subDir)

    elif peak_calling == 'Genrich':
        PEAK_dir = '%s/Peaks_Genrich' % root_dir
        QTL_dir = '%s/QTLs_Genrich/%s' % (root_dir, GT_subDir)

    elif peak_calling == 'Genrich_combined':
        PEAK_dir = '%s/Peaks_Genrich/combined' % root_dir
        QTL_dir = '%s/QTLs_Genrich_combined/%s' % (root_dir, GT_subDir)

    real_QTL_sig = compute_QTL_peakLevel(WINDOW, QTL_dir, saveSuffix = '_realGT_all')
    call_QTL_sig = compute_QTL_peakLevel(WINDOW, QTL_dir, saveSuffix = '_Called')
    call_QTL_imputed_sig = compute_QTL_peakLevel(WINDOW, '%s/Imputation' % QTL_dir, saveSuffix = '_Imputation')
    call_QTL_integrated_sig = compute_QTL_peakLevel(WINDOW, '%s/Integration' % QTL_dir, saveSuffix = '_Integration')
