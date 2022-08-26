[![CI](https://github.com/pontushojer/NAIBR/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/pontushojer/NAIBR/actions/workflows/ci.yml) [![codecov](https://codecov.io/gh/pontushojer/NAIBR/branch/main/graph/badge.svg?token=lHWwvvKaz1)](https://codecov.io/gh/pontushojer/NAIBR)
# NAIBR - Novel Adjacency Identification with Barcoded Reads

*NB! This is a re-factoring of the code from [raphael-group/NAIBR](https://github.com/raphael-group/NAIBR) to make the tool a bit more user-friendly and adds some additional features. This means that the resulting output might differ from what you would get from `raphael-group/NAIBR`.*

## Overview

NAIBR (Novel Adjacency Identification with Barcoded Reads) identifies novel adjacencies created by structural variation events such as deletions, duplications, inversions, and complex rearrangements using linked-read whole-genome sequencing data produced by 10X Genomics. Please refer to the [publication](https://doi.org/10.1093/bioinformatics/btx712) for details about the method.


NAIBR takes as in put a phased BAM, such as those produced by 10X Genomic's 
Long Ranger pipeline, and outputs a BEDPE file containing predicted novel 
adjacencies and a likelihood score for each adjacency.

## Installing NAIBR
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

* bam_file: Input BAM file < required >
* min_mapq: Minimum mapping quality for a read to be included in analysis (default: 40)
* outdir: Output directory (default: . )
* d: The maximum distance between reads in a linked-read
* blacklist: tap separated list of regions to be excluded from analysis (default: None)
* candidates: List in BEDPE format of novel adjacencies to be scored by NAIBR. This will override automatic detection of candidate novel adjacencies. 
* threads: Number of threads (default: 1)
* min_sv: Minimum size of a structural variant to be detected (default: lmax, the 95th percentile of the paired-end read insert size distribution)
* k: minimum number of barcode overlaps supporting a candidate NA (default = 3)

## Output

Scored novel adjacencies are outputted in three different file formats to the output directory

### `NAIBR.bedpe` 
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

### `NAIBR.reformat.bedpe`
Reformatted BEDPE that follows the [10x Genomics BEDPE format](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bedpe) suitable for IGV visualization


### `NAIBR.vcf`

NAs translated to the VCF format. Useful for operation using e.g. [truvari](https://github.com/ACEnglish/truvari).

## Example
Example files are provided in the 'example' directory. Running

```
naibr example/example.config
```

will produce the file 'example/NAIBR_SVs.bedpe'.

## Converting NAIBR BEDPE to VCF

The repository includes a utility script `scripts/bedpe_to_vcf.py` for converting from NAIBR's BEDPE format to a VCF. For detail on how to use it run:

```
python scripts/bedpe_to_vcf.py --help
```

### Citing NAIBR
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


