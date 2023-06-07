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
