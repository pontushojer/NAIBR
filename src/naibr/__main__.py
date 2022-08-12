"""
NAIBR - Novel Adjacency Identification with Barcoded Reads

Usage: naibr input.config

Running requires a config file with "=" separated parameters and values.

required:
    bam_file - str - Path to BAM file with phased linked reads, i.e. has BX and HP tags.

options:
    outdir - str - Path to output directory. Default: CWD
    blacklist - str - Path to BED file with regions to ignore
    candidates - str - Path to BEDPE with candidiate SVs to evaluate
    d - int - Maximum distance between read-pairs in a linked-read. Default: 10,000
    min_mapq - int - Minimum mapping quality to evaluate reads. Default: 40
    k - int - Minimum number of barcode overlaps supporting a candidate NA. Default: 3
    min_sv - int - Minimum size of structural variant. Default: 1000
    sd_mult - int - Stddev multiplier for max/min insert size (LMAX/LMIN) calculcation. Default: 2
    min_len - int - Minimum length of linked read fragment. Default: 2*LMAX
    max_len - int - TO_BE_ADDED. Default: 200,000
    min_reads - int - Minimum reads in linked read fragment. Default: 2
    min_discs - int - Minimum number of discordant reads. Default: 2
    threads - int - Threads to use. Default: 1
    DEBUG - bool - Run in debug mode. Default: False

About:

    NAIBR identifies novel adjacencies created by structural variation events such as
    deletions, duplications, inversions, and complex rearrangements using linked-read
    whole-genome sequencing data as produced by 10X Genomics. Please refer to the
    publication for details about the method.

Citation:

    Elyanow R, Wu HT, Raphael BJ. Identifying structural variants using linked-read
    sequencing data. Bioinformatics. 2018 Jan 15;34(2):353-360.
    doi: 10.1093/bioinformatics/btx712
"""
import sys
import time
import pysam
import collections
import numpy as np

if len(sys.argv) < 2 or sys.argv[1] in {"help", "-h", "--help"}:
    sys.exit(__doc__)

from .get_reads import parse_candidate_region, get_distributions, parse_chromosome, get_candidates
from .global_vars import configs
from .utils import flatten, parallel_execute, is_proper_chrom, write_novel_adjacencies
from .rank import predict_novel_adjacencies


def run_naibr_candidate(cand):
    """
    use user input candidate novel adjacencies
    """
    scores = 0
    reads_by_barcode, barcodes_by_pos, discs_by_barcode, cands, coverage = parse_candidate_region(cand)

    # Filter out barcodes with too few reads
    if configs.MIN_READS > 1:
        reads_by_barcode = {k: v for k, v in reads_by_barcode.items() if len(v) >= configs.MIN_READS}

    p_len, p_rate, overlap = get_distributions(reads_by_barcode)
    if p_len is None:
        return scores
    scores = predict_novel_adjacencies(
        reads_by_barcode, barcodes_by_pos, discs_by_barcode, cands, p_len, p_rate, coverage, False
    )
    return scores


def run_naibr(chrom):
    """
    automatically identify candidate novel adjacencies
    """
    (
        reads_by_barcode,
        barcodes_by_pos,
        discs_by_barcode,
        discs,
        interchrom_discs,
        coverage,
    ) = parse_chromosome(chrom)

    # Filter out barcodes with too few reads
    if configs.MIN_READS > 1:
        reads_by_barcode = {k: v for k, v in reads_by_barcode.items() if len(v) >= configs.MIN_READS}

    novel_adjacencies = []
    if len(reads_by_barcode) > 0:
        t_get_candidates = time.time()
        cands, p_len, p_rate = get_candidates(discs, reads_by_barcode)
        print(f"Got candidate noval adjacencies from data: {time.time() - t_get_candidates:.4f} s")
        if cands is None:
            print("No candidates from %s" % chrom)
        else:
            print("ranking %i candidates from %s" % (len(cands), chrom))
            novel_adjacencies = predict_novel_adjacencies(
                reads_by_barcode,
                barcodes_by_pos,
                discs_by_barcode,
                cands,
                p_len,
                p_rate,
                coverage,
                False,
            )
    else:
        print("No candidates from %s" % chrom)
    return reads_by_barcode, barcodes_by_pos, discs_by_barcode, interchrom_discs, coverage, novel_adjacencies


def main():
    starttime = time.time()

    print("========= NAIBR =========")
    print("-------- CONFIGS --------")
    print("FILES")
    print(f"  bam_file = {configs.BAM_FILE}")
    print(f"  candidates = {configs.CANDIDATES}")
    print(f"  blacklist = {configs.BLACKLIST_FILE}")
    print(f"  outdir = {configs.DIR}")
    print("PARAMETERS")
    print(f"  d =         {configs.MAX_LINKED_DIST:>9,}")
    print(f"  min_mapq =  {configs.MIN_MAPQ:>9}")
    print(f"  k =         {configs.MIN_BC_OVERLAP:>9}")
    print(f"  min_sv =    {configs.MIN_SV:>9,}")
    print(f"  sd_mult =   {configs.SD_MULT:>9}")
    print(f"  min_len =   {configs.MIN_LEN:>9,}")
    print(f"  max_len =   {configs.MAX_LEN:>9,}")
    print(f"  min_reads = {configs.MIN_READS:>9}")
    print(f"  min_discs = {configs.MIN_DISCS:>9}")
    print(f"  threads =   {configs.NUM_THREADS:>9}")
    print(f"  DEBUG =     {str(configs.DEBUG).rjust(9)}")
    print("-------------------------")

    if configs.CANDIDATES:
        print("Using user defined candidates")
        cands = []
        with open(configs.CANDIDATES) as f:
            for line in f:
                els = line.strip().split("\t")
                if len(els) > 4:
                    cands.append([els[0], int(els[1]), els[3], int(els[4]), els[-1]])
        novel_adjacencies = flatten(parallel_execute(run_naibr_candidate, cands))
        write_novel_adjacencies(novel_adjacencies)

    else:
        chroms = []
        with pysam.AlignmentFile(configs.BAM_FILE, "rb") as reads:
            chroms = reads.references
            chroms = [x for x in chroms if is_proper_chrom(x)]

        data = parallel_execute(run_naibr, chroms)
        reads_by_barcode = collections.defaultdict(list)
        barcodes_by_pos = collections.defaultdict(list)
        discs_by_barcode = collections.defaultdict(list)
        discs = collections.defaultdict(list)
        coverage = []
        novel_adjacencies = []
        for (
            reads_by_barcode_chrom,
            barcodes_by_pos_chrom,
            discs_by_barcode_chrom,
            discs_chrom,
            cov_chrom,
            nas_chrom,
        ) in data:
            if nas_chrom:
                reads_by_barcode.update(reads_by_barcode_chrom)
                barcodes_by_pos.update(barcodes_by_pos_chrom)
                discs_by_barcode.update(discs_by_barcode_chrom)
                discs.update(discs_chrom)
                coverage.append(cov_chrom)
                novel_adjacencies += nas_chrom
        cands, p_len, p_rate = get_candidates(discs, reads_by_barcode)
        if cands is not None:
            print("ranking %i interchromosomal candidates" % len(cands))
            novel_adjacencies += predict_novel_adjacencies(
                reads_by_barcode,
                barcodes_by_pos,
                discs_by_barcode,
                cands,
                p_len,
                p_rate,
                np.mean(coverage),
                True,
            )
        else:
            print("No interchromosomal candidates")
        write_novel_adjacencies(novel_adjacencies)

    print("Finished in", (time.time() - starttime) / 60.0, "minutes")
    return


if __name__ == "__main__":
    main()
