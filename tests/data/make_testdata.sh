#!/bin/bash
#
# Generate test data for the NAIBR test suite, `samtools` required to run
#

set -euox pipefail

# TRANSLOCATION ON CHROMOSOME 12 and chr20 (HCC1954T)
# hg19: chr12	40393752	40393756	-> chr20	53974890	53974896
# Source: https://doi.org/10.1093/nar/gkr221
bam="https://s3-us-west-2.amazonaws.com/10x.files/samples/genome/HCC1954T_WGS_210/HCC1954T_WGS_210_phased_possorted_bam.bam"
region1="chr12:40,343,752-40,433,756"
region2="chr20:53,934,890-54,024,896"

samtools view -bh "${bam}" "${region1}" -o HCC1954T_10xGenomics_chr_12_TRA.tmp.bam
samtools view -bh "${bam}" "${region2}" -o HCC1954T_10xGenomics_chr_20_TRA.tmp.bam
samtools cat -o HCC1954T_10xGenomics_chr_12_20_TRA.bam HCC1954T_10xGenomics_chr_12_TRA.tmp.bam HCC1954T_10xGenomics_chr_20_TRA.tmp.bam
samtools index HCC1954T_10xGenomics_chr_12_20_TRA.bam
rm HCC1954T_10xGenomics_chr_20_TRA.tmp.bam HCC1954T_10xGenomics_chr_12_TRA.tmp.bam

# INVERSION ON CHROMOSOME 21 (HCC1954T)
# This is the same as used in the https://github.com/raphael-group/NAIBR/ example
# hg19: chr21   18759465    18906782
bam="https://s3-us-west-2.amazonaws.com/10x.files/samples/genome/HCC1954T_WGS_210/HCC1954T_WGS_210_phased_possorted_bam.bam"
region="chr21:18,709,465-18,956,782"
samtools view "${bam}" "${region}" -o HCC1954T_10xGenomics_chr21_INV.bam
samtools index HCC1954T_10xGenomics_chr21_INV.bam

# DELETION ON CHROMOSOME 1 (NA24385)
# GRCh38: chr1  115686862    115690219
# source: CMRG benchmark,  https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/
bam="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh38/NA24385.GRCh38.phased_possorted_bam.bam"
region="chr1:115,620,601-115,756,835"
samtools view "${bam}" "${region}" -o NA24385_10xGenomics_chr1_DEL.bam
samtools index NA24385_10xGenomics_chr1_DEL.bam

# DUPLICATION ON CHROMOSOME 3 (HCC1954T)
# GRCh38: chr3  101569210   101663759
# source: LongRanger calls
bam="https://s3-us-west-2.amazonaws.com/10x.files/samples/genome/HCC1954T_WGS_210/HCC1954T_WGS_210_phased_possorted_bam.bam"
region="chr3:101,508,502-101,725,227"
samtools view "${bam}" "${region}" -o HCC1954T_10xGenomics_chr3_DUP.bam
samtools index HCC1954T_10xGenomics_chr3_DUP.bam


# Subsample BAM to contain only a single barcode
samtools view -d BX:AACGGTTCACGAAACG-1 tests/data/HCC1954T_10xGenomics_chr21_INV.bam -bh -o tests/data/HCC1954T_10xGenomics_chr21_INV.one_barcode.bam
samtools index tests/data/HCC1954T_10xGenomics_chr21_INV.one_barcode.bam
