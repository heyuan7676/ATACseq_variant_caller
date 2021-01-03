import sys
import os
import pandas as pd


if __name__ == "__main__":
    ROOT_DIR = sys.argv[1]  # /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/
    CHROMOSOME = int(sys.argv[2])
    peak_calling = sys.argv[3]
    WINDOW = int(sys.argv[4])
    minDP = int(sys.argv[5])

    alignment_dir = 'alignment_bowtie' # alignment_subsample_0.5
    GT_subDir = 'minDP%d' % minDP

    root_dir = '%s/%s' % (ROOT_DIR, alignment_dir)
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'

    if peak_calling == 'macs2':
    	PEAK_dir = '%s/Peaks_MACS2' % root_dir
        QTL_dir = '%s/QTLs_MACS2/%s' % (root_dir, GT_subDir)

    if peak_calling == 'macs2_combined':
        PEAK_dir = '%s/Peak_version1/Peaks/combined' % root_dir
        QTL_dir = '%s/QTLs_MACS2_combined/%s' % (root_dir, GT_subDir)

    elif peak_calling == 'Genrich':
        PEAK_dir = '%s/Peaks_Genrich' % root_dir
        QTL_dir = '%s/QTLs_Genrich/%s' % (root_dir, GT_subDir)

    elif peak_calling == 'Genrich_combined':
        PEAK_dir = '%s/Peaks_Genrich/combined' % root_dir
        QTL_dir = '%s/QTLs_Genrich_combined/%s' % (root_dir, GT_subDir)

    if not os.path.exists(QTL_dir):
        os.makedirs(QTL_dir)

    QTL_dir_imputed = '%s/Imputation' % QTL_dir
    if not os.path.exists(QTL_dir_imputed):
	os.makedirs(QTL_dir_imputed)

    QTL_dir_integration = '%s/Integration' % QTL_dir
    if not os.path.exists(QTL_dir_integration):
        os.makedirs(QTL_dir_integration)

    fn = os.path.join(QTL_dir_integration, 'CHR%d_caQTLs_WINDOW_%skb%s.txt' % (CHROMOSOME, str(WINDOW/1000.0), '_Integration'))
    flag = 0
    try:
        dat = pd.read_csv(fn, nrows=10, sep='\t')
        flag = 1
    except:
        print('Call ca-QTLs for chromosome %d with window = %skb, %s - Integration' % (CHROMOSOME, str(WINDOW/1000.0), peak_calling))


    fn = os.path.join(QTL_dir_imputed, 'CHR%d_caQTLs_WINDOW_%skb%s.txt' % (CHROMOSOME, str(WINDOW/1000.0), '_Imputation'))
    flag = 0
    try:
        dat = pd.read_csv(fn, nrows=10, sep='\t')
        flag = 1
    except:
        print('Call ca-QTLs for chromosome %d with window = %skb, %s - Imputation' % (CHROMOSOME, str(WINDOW/1000.0), peak_calling))
    

    fn = os.path.join(QTL_dir, 'CHR%d_caQTLs_WINDOW_%skb%s.txt' % (CHROMOSOME, str(WINDOW/1000.0), '_Called'))
    flag = 0
    try:
        dat = pd.read_csv(fn, nrows=10, sep='\t')
        flag = 1
    except:
        print('Call ca-QTLs for chromosome %d with window = %skb, %s - Called' % (CHROMOSOME, str(WINDOW/1000.0), peak_calling))


    fn = os.path.join(QTL_dir, 'CHR%d_caQTLs_WINDOW_%skb%s.txt' % (CHROMOSOME, str(WINDOW/1000.0), '_realGT_all'))
    flag = 0
    try:
        dat = pd.read_csv(fn, nrows=10, sep='\t')
        flag = 1
    except:
        print('Call ca-QTLs for chromosome %d with window = %skb, %s - WGS' % (CHROMOSOME, str(WINDOW/1000.0), peak_calling))

