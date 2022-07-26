### GLOBALS ####
from __future__ import print_function
import os
import sys
import collections
import pysam
import numpy as np


def parse_blacklist(fname):
    blacklist = collections.defaultdict(list)
    with open(fname) as f:
        for line in f:
            fields = line.strip().split("\t")
            assert len(fields) == 3
            chrom, start, end = fields
            chrom = chrom.strip("chr")
            blacklist[chrom].append([int(start), int(end)])
    return blacklist


def parse_configs(fname):
    configs = {}
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue

            name, val = line.split("=")
            configs[name] = val
    return configs


constants = parse_configs(sys.argv[1])


if "min_mapq" in constants:
    MIN_MAPQ = int(constants["min_mapq"])
else:
    MIN_MAPQ = 40

if "k" in constants:
    MIN_BC_OVERLAP = int(constants["k"])
else:
    MIN_BC_OVERLAP = 3

if "bam_file" in constants and os.path.exists(constants["bam_file"]):
    BAM_FILE = constants["bam_file"]
else:
    sys.exit("Missing path to BAM file in config file. Please include this with the "
             "pattern 'bam_file=/path/to/file.bam'")

if "DEBUG" in constants:
    DEBUG = bool(constants["DEBUG"])
else:
    DEBUG = False

if "outdir" in constants:
    DIR = constants["outdir"]
    if not os.path.exists(DIR):
        os.makedirs(DIR)
else:
    DIR = ".."

if "d" in constants:
    MAX_LINKED_DIST = int(constants["d"])
else:
    MAX_LINKED_DIST = 10_000

if "threads" in constants:
    NUM_THREADS = int(constants["threads"])
else:
    NUM_THREADS = 1

if "sd_mult" in constants:
    SD_MULT = int(constants["sd_mult"])
else:
    SD_MULT = 2


def parse_read_pairs(iterator):
    mates = {}
    for read in iterator:
        if read.is_unmapped or read.mate_is_unmapped:
            continue

        if read.query_name in mates:
            yield read, mates.pop(read.query_name)
        else:
            mates[read.query_name] = read


def estimate_lmin_lmax():
    reads = pysam.AlignmentFile(BAM_FILE, "rb")
    pair_spans = []
    reads_lengths = []
    num = 0
    for i, chrm in enumerate(reads.references):
        for read, mate in parse_read_pairs(reads.fetch(chrm)):
            num += 1
            if num > 1_000_000:
                break

            start = min(read.reference_start, mate.reference_start)
            end = max(read.reference_end, mate.reference_end)
            dist = end-start

            if 0 < abs(dist) < 2000:
                reads_lengths.extend([read.query_length, mate.query_length])
                pair_spans.append(dist)

    mean_dist = np.mean(pair_spans)
    std_dist = np.std(pair_spans)
    lmin = max(int(mean_dist - std_dist * SD_MULT), -int(np.mean(reads_lengths)))
    lmax = max(int(mean_dist + std_dist * SD_MULT), 100)

    if DEBUG:
        print(f"lmax = {lmax}, lmin = {lmin}")
    return lmin, lmax


LMIN, LMAX = estimate_lmin_lmax()


if "min_sv" in constants:
    MIN_SV = int(constants["min_sv"])
else:
    MIN_SV = 2 * LMAX

if "min_len" in constants:
    MIN_LEN = int(constants["min_len"])
else:
    MIN_LEN = 2 * LMAX

if "max_len" in constants:
    MAX_LEN = int(constants["max_len"])
else:
    MAX_LEN = 200000

if "min_reads" in constants:
    MIN_READS = int(constants["min_reads"])
else:
    MIN_READS = 2

if "min_discs" in constants:
    MIN_DISCS = int(constants["min_discs"])
else:
    MIN_DISCS = 2

if "candidates" in constants and os.path.exists(constants["candidates"]):
    CANDIDATES = constants["candidates"]
else:
    CANDIDATES = ""

if "blacklist" in constants and os.path.exists(constants["blacklist"]):
    BLACKLIST_FILE = constants["blacklist"]
    BLACKLIST = parse_blacklist(constants["blacklist"])
else:
    BLACKLIST_FILE = ""
    BLACKLIST = None
