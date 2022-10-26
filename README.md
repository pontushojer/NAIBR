[![CI](https://github.com/pontushojer/NAIBR/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/pontushojer/NAIBR/actions/workflows/ci.yml) [![codecov](https://codecov.io/gh/pontushojer/NAIBR/branch/main/graph/badge.svg?token=lHWwvvKaz1)](https://codecov.io/gh/pontushojer/NAIBR) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
# NAIBR - Novel Adjacency Identification with Barcoded Reads

> **Note** This is a re-factoring of the code from [raphael-group/NAIBR](https://github.com/raphael-group/NAIBR) to make the tool a bit more user-friendly and adds some additional features. This means that the resulting output might differ from what you would get from `raphael-group/NAIBR`. See [Benchmarks](#Benchmarks) for a comparison. 

## Contents

- [Overview](#Overview)
- [Install](#Install)
- [Running NAIBR](#Running-NAIBR)
- [Converting BEDPE to VCF](#Converting-BEDPE-to-VCF)
- [Benchmarks](#Benchmarks)
- [Citing](#Citing)

## Overview

NAIBR (Novel Adjacency Identification with Barcoded Reads) identifies novel adjacencies (NAs) created by structural variation events such as deletions, duplications, inversions, and complex rearrangements using linked-read whole-genome sequencing data produced by 10X Genomics. Please refer to the [publication](https://doi.org/10.1093/bioinformatics/btx712) for details about the method.


NAIBR takes as in put a phased BAM, such as those produced by 10X Genomic's 
Long Ranger pipeline, and outputs a BEDPE file containing predicted novel 
adjacencies and a likelihood score for each adjacency.

## Install
```
git clone https://github.com/pontushojer/NAIBR.git
cd NAIBR
pip install .
```

## Running NAIBR

NAIBR can be run using the following command:

```
naibr <configfile>
```

A template config file can be found in example/example.config. The following parameters can be set in the config file:

* `bam_file`: Input phased BAM file with `BX` and `HP` tagged reads (**REQUIRED**)
* `min_mapq`: Minimum mapping quality for a read to be included in analysis (default: 40)
* `outdir`: Output directory (default: `.` i.e. current workding directory)
* `d`: The maximum distance in basepairs between reads in a linked-read (default: 10000)
* `blacklist`: BED-file with regions to be excluded from analysis
* `candidates`: BEDPE-file with novel adjacencies to be scored by NAIBR. This will override automatic detection of candidate novel adjacencies 
* `threads`: Number of threads (default: 1)
* `min_sv`: Minimum size of a structural variant to be detected (default: `lmax`, i.e. the 95th percentile of the paired-end read insert size distribution)
* `k`: minimum number of barcode overlaps supporting a candidate NA (default = 3)


### Run example
Example files are provided in the 'example' directory. Running

```
naibr example/example.config
```

will produce the files `example/NAIBR_SVs.bedpe`, `example/NAIBR_SVs.reformat.bedpe` and `example/NAIBR_SVs.vcf`.

### Output

Scored novel adjacencies are outputted in three different file formats to the output directory

#### `NAIBR.bedpe` 
BEDPE-like file in same format as outputted in [raphael-group/NAIBR](https://github.com/raphael-group/NAIBR). Note that this file does not follow any actual BEDPE format. The file has a header with column descriptions.

**Columns:**
1. Chr1: Chromosome of 1st breakpoint
2. Break1: Position of 1st breakpoint
3. Chr2: Chromosome of 2nd breakpoint
4. Break2: Position of 2nd breakpoint
5. Split molecules: Number of split molecule supporting NA
6. Discordant reads: Number of discordant reads supporting NA 
7. Orientation: Orientation of NA. "+-" ~ DEL, "++"/"--" ~ INV, "-+" ~ "DUP"
8. Haplotype
9. Score: log-likelihood score for NA
10. Pass filter: "PASS" if passes filter threshold, otherwice "FAIL"

#### `NAIBR.reformat.bedpe`
Reformatted BEDPE that follows the [10x Genomics BEDPE format](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bedpe) suitable for IGV visualization


#### `NAIBR.vcf`

NAs translated to the VCF format. Useful for operation using e.g. [truvari](https://github.com/ACEnglish/truvari).

## Converting BEDPE to VCF

The repository includes a utility script `scripts/bedpe_to_vcf.py` for converting from NAIBR's BEDPE format to a VCF. For detail on how to use it run:

```
python scripts/bedpe_to_vcf.py --help
```

## Benchmarks
Benchmarking this fork (labeled `main_<commit>`) versus [raphael-group/NAIBR](https://github.com/raphael-group/NAIBR) (labeled `master`). Master includes https://github.com/raphael-group/NAIBR/pull/22 merged as this removes a lot of FP calls introduced by a code error.

### HG002 
Run NAIBR on 10x Genomics data for idividual HG002 (NA24385). Dataset from GIAB ([source](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh37/)) has 51.7x coverage and was run through the Long Ranger pipeline using reference genome GRCh37. 

<details>
  <summary>NAIBR configs used</summary>
  <code>
  min_mapq=40
  d=10000
  min_sv=1000
  threads=10
  k=3
  min_reads=2
  min_discs=2
</code>
</details>

Compared deletions (DELs) to [GIAB benchmark](https://doi.org/10.1038%2Fs41587-020-0538-8) set `Tier1` (v0.6) using [truvari](https://github.com/ACEnglish/truvari) (v3.1.0). Compared both all entries (named `*_all`) and only those that passed filters (labeled `*_pass`).

Command for `*_all`:
```
truvari bench -b HG002_SVs_Tier1_v0.6.withchr.DELs.vcf.gz -c master/NAIBR_SVs.DELs.vcf.gz -o truvari_bench/master_all --includebed HG002_SVs_Tier1_v0.6.withchr.noY.bed --pctsim 0.5 --pctsize 0.5 --reference /Users/pontus.hojer/analysis/references/hg19_longranger/genome.fa --sizemin 1000 --sizemax 20000000 --multimatch --giabreport
```

Command for `*_pass`:
```
truvari bench -b HG002_SVs_Tier1_v0.6.withchr.DELs.vcf.gz -c master/NAIBR_SVs.DELs.vcf.gz -o truvari_bench/master_pass --includebed HG002_SVs_Tier1_v0.6.withchr.noY.bed --pctsim 0.5 --pctsize 0.5 --reference /Users/pontus.hojer/analysis/references/hg19_longranger/genome.fa --sizemin 1000 --sizemax 20000000 --multimatch --giabreport --passonly
```

**Results**
Sample | Precision | Recall | F1 score | TP-base | TP-call | FP | FN
--- | --- | --- | --- | --- | --- | --- | ---
`main_cbbdf01_all`  | 67.0% | 59.6% | 63.1% | **480** | **478** | 235 | 325
`main_cbbdf01_pass` | **99.3%** | **77.4%** | **87.0%** | 411 | 411 | **3** | **120**
`master_all`        | 34.3% | 28.6% | 31.2% | 230 | 229 | 439 | 575
`master_pass`       | 96.3% | 14.9% | 25.8% | 79 | 79 | **3** | 452

This fork, i.e. `main_cbbdf01_pass`, had much higher F1 score of 87.0% than the master branch at 25.8%. Also the `main_cbbdf01_all` seems to include a lot of TP calls that are filtered out by NAIBR, possibly filtering should be re-optimized.

### Citing
Elyanow, Rebecca, Hsin-Ta Wu, and Benjamin J. Raphael. *Identifying 
structural variants using linked-read sequencing data.* Bioinformatics (2017).
```
@article{elyanow2017identifying,
  title={Identifying structural variants using linked-read sequencing data},
  author={Elyanow, Rebecca and Wu, Hsin-Ta and Raphael, Benjamin J},
  journal={Bioinformatics},
  year={2017}
}
```


