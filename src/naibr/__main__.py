"""
NAIBR - Novel Adjacency Identification with Barcoded Reads

Usage: naibr input.config

Running requires a config file with "=" separated parameters and values.

required:
    bam_file   - str  - Path to BAM file with phased linked reads, i.e. has BX and HP tags.

options:
    outdir     - str  - Path to output directory. Default: CWD
    blacklist  - str  - Path to BED file with regions to ignore
    candidates - str  - Path to BEDPE with candidiate SVs to evaluate
    d          - int  - Maximum distance between read-pairs in a linked-read. Default: 10,000
    min_mapq   - int  - Minimum mapping quality to evaluate reads. Default: 40
    k          - int  - Minimum number of barcode overlaps supporting a candidate NA. Default: 3
    min_sv     - int  - Minimum size of structural variant. Default: 1000
    sd_mult    - int  - Stddev multiplier for max/min insert size (LMAX/LMIN) calculcation. Default: 2
    min_len    - int  - Minimum length of linked read fragment. Default: 2*LMAX
    max_len    - int  - TO_BE_ADDED. Default: 200,000
    min_reads  - int  - Minimum reads in linked read fragment. Default: 2
    min_discs  - int  - Minimum number of discordant reads. Default: 2
    threads    - int  - Threads to use. Default: 1
    prefix     - str  - Prefix for output files. Default: NAIBR_SVs
    DEBUG      - bool - Run in debug mode. Default: False

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
import collections
import functools
import logging
import os
import sys
import time

import numpy as np
import pysam

from . import __version__
from .distributions import get_linked_reads
from .get_reads import get_candidates, get_interchrom_candidates, parse_candidate_region, parse_chromosome
from .global_vars import Configs
from .rank import predict_novel_adjacencies
from .utils import (
    ConditionalFormatter,
    flatten,
    input_candidates,
    is_proper_chrom,
    parallel_execute,
    UnionDict,
    write_novel_adjacencies,
)

logger = logging.getLogger(__name__)


def run_naibr_on_candidates(configs):
    """Run NAIBR on user input candidate novel adjacencies"""
    with open(configs.CANDIDATES) as f:
        cands = input_candidates(f, min_sv=configs.MIN_SV)

    if len(cands) == 0:
        logger.warning(f"No valid candidates in file {configs.CANDIDATES}.")
        return []

    configs.COMPRESSION_THREADS = 2 if configs.NUM_THREADS / len(cands) > 2 else 1
    run_with_configs = functools.partial(evaluate_candidate, configs=configs)
    novel_adjacencies = flatten(parallel_execute(run_with_configs, cands, threads=configs.NUM_THREADS))
    logger.info(f"Evaluation yielded {len(novel_adjacencies)} viable novel adjacencies")
    return novel_adjacencies


def evaluate_candidate(cand, configs):
    """
    use user input candidate novel adjacencies
    """
    (
        reads_by_barcode,
        barcodes_by_pos,
        discs_by_barcode,
        discs,
        cand,
        interchrom_discs,
        coverage,
    ) = parse_candidate_region(cand, configs)
    # Filter out barcodes with too few reads
    if configs.MIN_READS > 1:
        reads_by_barcode = {k: v for k, v in reads_by_barcode.items() if len(v) >= configs.MIN_READS}

    # See if there are other candidates in close proximity
    if not cand.is_interchromosomal():
        cand_name = f"{cand.chrm1}:{cand.break1}-{cand.break2}"
        logger.info(f"For {cand.svtype()} {cand_name} - found {len(discs):,} discordant positions")
        candidates, p_len, p_rate, _ = get_candidates(discs, reads_by_barcode, configs)
    else:
        cand_name = f"{cand.chrm1}:{cand.break1}-{cand.chrm2}:{cand.break2}"
        logger.info(f"For {cand.svtype()} {cand_name} - found {len(interchrom_discs):,} discordant positions")
        linkedreads_by_barcode = get_linked_reads(reads_by_barcode, configs)
        candidates, p_len, p_rate = get_interchrom_candidates(
            interchrom_discs, linkedreads_by_barcode, configs
        )

    if candidates:
        candidates = [c for c in candidates if max(c.distance(cand)) < configs.LMAX]
        candidates.append(cand)
    else:
        candidates = [cand]

    if p_len is None:
        return []

    logger.info(f"ranking {len(candidates)} candidates for {cand_name}")
    if not cand.is_interchromosomal():
        scores = predict_novel_adjacencies(
            reads_by_barcode,
            barcodes_by_pos,
            discs_by_barcode,
            candidates,
            p_len,
            p_rate,
            coverage,
            False,
            configs,
        )
    else:
        scores = predict_novel_adjacencies(
            linkedreads_by_barcode,
            barcodes_by_pos,
            discs_by_barcode,
            candidates,
            p_len,
            p_rate,
            coverage,
            True,
            configs,
        )
    return scores


def run_naibr_on_chromosomes(chromosomes, configs):
    # Setup
    configs.COMPRESSION_THREADS = 2 if configs.NUM_THREADS / len(chromosomes) > 2 else 1
    run_with_configs = functools.partial(run_naibr_on_chromosome, configs=configs)

    # Run
    chroms_data = parallel_execute(run_with_configs, chromosomes, threads=configs.NUM_THREADS)

    # Collect NAs and data from chromosomes
    linkedreads_by_barcode = UnionDict(list)
    barcodes_by_pos = UnionDict(set)
    discs_by_barcode = UnionDict(list)
    interchrom_discs = UnionDict(list)
    coverage = []
    novel_adjacencies = []
    for data in chroms_data:
        (
            linkedreads_by_barcode_chrom,
            barcodes_by_pos_chrom,
            discs_by_barcode_chrom,
            interchrom_discs_chrom,
            cov_chrom,
            nas_chrom,
        ) = data
        linkedreads_by_barcode.combine(linkedreads_by_barcode_chrom)
        barcodes_by_pos.combine(barcodes_by_pos_chrom)
        discs_by_barcode.combine(discs_by_barcode_chrom)
        interchrom_discs.combine(interchrom_discs_chrom)
        coverage.append(cov_chrom)
        novel_adjacencies += nas_chrom

    if len(chromosomes) == 1:
        return novel_adjacencies

    # Get interchromosomal candidates
    cands, p_len, p_rate = get_interchrom_candidates(interchrom_discs, linkedreads_by_barcode, configs)

    if not cands:
        logger.info("No interchromosomal candidates")
        return novel_adjacencies

    logger.info("ranking %i interchromosomal candidates" % len(cands))
    interchrom_novel_adjacencies = predict_novel_adjacencies(
        linkedreads_by_barcode,
        barcodes_by_pos,
        discs_by_barcode,
        cands,
        p_len,
        p_rate,
        np.mean(coverage),
        True,
        configs,
    )
    if interchrom_novel_adjacencies:
        n_pass = sum([na.pass_threshold for na in interchrom_novel_adjacencies])
        n_total = len(interchrom_novel_adjacencies)
        logger.info(
            f"Found {n_total:,} interchromosomal SVs of which {n_pass:,} "
            f"({n_pass / n_total:.1%}) passed the threshold"
        )
        novel_adjacencies += interchrom_novel_adjacencies

    return novel_adjacencies


def run_naibr_on_chromosome(chrom, configs):
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
    ) = parse_chromosome(chrom, configs)

    # Filter out barcodes with too few reads
    if configs.MIN_READS > 1:
        reads_by_barcode = {k: v for k, v in reads_by_barcode.items() if len(v) >= configs.MIN_READS}

    novel_adjacencies = []
    linkedreads_by_barcode = collections.defaultdict(list)
    if len(reads_by_barcode) > 0:
        t_get_candidates = time.time()
        logger.info(f"For chromsome {chrom} - found {len(discs):,} positions with discordant reads.")
        cands, p_len, p_rate, linkedreads_by_barcode = get_candidates(discs, reads_by_barcode, configs)
        logger.info(f"Got candidate noval adjacencies from data: {time.time() - t_get_candidates:.4f} s")
        if not cands:
            logger.info("No candidates from %s" % chrom)
        else:
            logger.info("Ranking %i candidates from %s" % (len(cands), chrom))
            novel_adjacencies = predict_novel_adjacencies(
                reads_by_barcode,
                barcodes_by_pos,
                discs_by_barcode,
                cands,
                p_len,
                p_rate,
                coverage,
                False,
                configs,
            )
            if novel_adjacencies:
                n_pass = sum([na.pass_threshold for na in novel_adjacencies])
                n_total = len(novel_adjacencies)
                logger.info(
                    f"Found {n_total:,} SVs on {chrom} of which {n_pass:,} "
                    f"({n_pass / n_total:.1%}) passed the threshold"
                )

    else:
        logger.info("No candidates from %s" % chrom)
    return (
        linkedreads_by_barcode,
        barcodes_by_pos,
        discs_by_barcode,
        interchrom_discs,
        coverage,
        novel_adjacencies,
    )


def chromosome_has_barcoded_reads(chromosome, bam_file, max_iter=100_000):
    """True if found BX tagged read on chromsome within the first `max_iter` reads."""
    with pysam.AlignmentFile(bam_file, "rb") as reads:
        for nr, read in enumerate(reads.fetch(chromosome)):
            if read.has_tag("BX"):
                return True

            if nr > max_iter:
                return False
    return False


def get_chromosomes_with_reads(bam_file):
    """Returns a list of chromosomes/contigs which has reads containing BX tags"""
    with pysam.AlignmentFile(bam_file, "rb") as reads:
        all_chroms = reads.references
        all_chroms = [x for x in all_chroms if is_proper_chrom(x)]

    chroms = [chrom for chrom in all_chroms if chromosome_has_barcoded_reads(chrom, bam_file)]
    return chroms


def parse_args(args):
    """Parse command line arguments"""
    if len(args) != 1 or args[0] in {"help", "-h", "--help"}:
        print(__doc__)
        sys.exit(0)

    if args[0] in {"version", "-v", "--version"}:
        print(f"naibr {__version__}")
        sys.exit(0)

    if not os.path.exists(args[0]):
        logger.error(f"Configs file '{args[0]}' does not exist")
        sys.exit(1)

    with open(args[0]) as f:
        file_configs = Configs.from_file(f)

    return file_configs


def run(configs):
    starttime = time.time()

    logger.info("========= NAIBR =========", extra={"nofmt": True})
    logger.info(f"version: {__version__}", extra={"nofmt": True})
    logger.info("-------- CONFIGS --------", extra={"nofmt": True})
    logger.info("FILES", extra={"nofmt": True})
    logger.info(f"  bam_file = {configs.BAM_FILE}", extra={"nofmt": True})
    logger.info(f"  candidates = {configs.CANDIDATES}", extra={"nofmt": True})
    logger.info(f"  blacklist = {configs.BLACKLIST_FILE}", extra={"nofmt": True})
    logger.info(f"  outdir = {configs.DIR}", extra={"nofmt": True})
    logger.info(f"  prefix = {configs.PREFIX}", extra={"nofmt": True})
    logger.info("PARAMETERS", extra={"nofmt": True})
    logger.info(f"  d =         {configs.MAX_LINKED_DIST:>9,}", extra={"nofmt": True})
    logger.info(f"  min_mapq =  {configs.MIN_MAPQ:>9}", extra={"nofmt": True})
    logger.info(f"  k =         {configs.MIN_BC_OVERLAP:>9}", extra={"nofmt": True})
    logger.info(f"  min_sv =    {configs.MIN_SV:>9,}", extra={"nofmt": True})
    logger.info(f"  sd_mult =   {configs.SD_MULT:>9}", extra={"nofmt": True})
    logger.info(f"  min_len =   {configs.MIN_LEN:>9,}", extra={"nofmt": True})
    logger.info(f"  max_len =   {configs.MAX_LEN:>9,}", extra={"nofmt": True})
    logger.info(f"  min_reads = {configs.MIN_READS:>9}", extra={"nofmt": True})
    logger.info(f"  min_discs = {configs.MIN_DISCS:>9}", extra={"nofmt": True})
    logger.info(f"  threads =   {configs.NUM_THREADS:>9}", extra={"nofmt": True})
    logger.info(f"  DEBUG =     {str(configs.DEBUG).rjust(9)}", extra={"nofmt": True})
    logger.info("-------------------------", extra={"nofmt": True})

    novel_adjacencies = []
    if configs.CANDIDATES:
        logger.info("Using user defined candidates")
        novel_adjacencies = run_naibr_on_candidates(configs)
    else:
        chromosomes = get_chromosomes_with_reads(bam_file=configs.BAM_FILE)
        if len(chromosomes) == 0:
            logger.error("BAM does not contain BX tagged reads.")
            sys.exit(1)

        logger.info(
            f"Found {len(chromosomes)} chromosome{'s' if len(chromosomes) > 1 else ''} with data in BAM"
        )

        novel_adjacencies = run_naibr_on_chromosomes(chromosomes, configs)

    write_novel_adjacencies(
        novel_adjacencies, directory=configs.DIR, bam_file=configs.BAM_FILE, prefix=configs.PREFIX
    )

    logger.info(f"Finished in {(time.time() - starttime) / 60.0} minutes")
    return 0


def main(args=None):
    if args is None:
        args = sys.argv[1:] if len(sys.argv) > 1 else []

    file_configs = parse_args(args)

    # Setup logging
    root = logging.getLogger()
    root.setLevel(logging.INFO)

    # Special formatter that ca be turned off
    formatter = ConditionalFormatter("%(asctime)s - %(levelname)s: %(message)s")

    # This prints to stdout
    streamhandler = logging.StreamHandler(stream=sys.stdout)
    streamhandler.setFormatter(formatter)
    root.addHandler(streamhandler)

    # Any existing log file is overwritten
    log_file = file_configs.DIR / f"{file_configs.PREFIX}.log"
    if log_file.exists():
        os.remove(log_file)

    # This prints to a log file in the output directory
    filehandler = logging.FileHandler(log_file)
    filehandler.setFormatter(formatter)
    root.addHandler(filehandler)

    if file_configs.DEBUG:
        root.setLevel(logging.DEBUG)

    return run(file_configs)


if __name__ == "__main__":
    sys.exit(main())
