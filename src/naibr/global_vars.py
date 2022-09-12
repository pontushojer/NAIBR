import collections
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pysam


def parse_blacklist(fname):
    blacklist = collections.defaultdict(list)
    with open(fname) as f:
        for line in f:
            fields = line.strip().split("\t")
            assert len(fields) == 3
            chrom, start, end = fields
            blacklist[chrom].append([int(start), int(end)])
    return blacklist


def parse_read_pairs(iterator):
    mates = {}
    for read in iterator:
        if read.is_unmapped or read.mate_is_unmapped:
            continue

        if read.query_name in mates:
            yield read, mates.pop(read.query_name)
        else:
            mates[read.query_name] = read


def estimate_lmin_lmax(bam_file, sd_mult, debug=False):
    """
    Estmate the upper and lower bounds for insert sizes by looing at up to
    the first million read pairs.
    """
    reads = pysam.AlignmentFile(bam_file, "rb")
    pair_spans = []
    reads_lengths = []
    num = 0
    for chrm in reads.references:
        for read, mate in parse_read_pairs(reads.fetch(chrm)):
            num += 2
            if num > 1_000_000:
                break

            start = min(read.reference_start, mate.reference_start)
            end = max(read.reference_end, mate.reference_end)
            dist = end - start

            if abs(dist) < 2000:
                reads_lengths.extend([read.query_length, mate.query_length])
                pair_spans.append(dist)

    pair_spans = np.array(pair_spans)
    mean_dist = pair_spans.mean()
    std_dist = pair_spans.std()

    lmin = max(int(mean_dist - std_dist * sd_mult), int(np.mean(reads_lengths)))
    lmax = max(int(mean_dist + std_dist * sd_mult), 100)

    if debug:
        # Estimate how large a percentage of reads have span within [lmin, lmax].
        reads_included = len(pair_spans[(pair_spans < lmax) & (pair_spans > lmin)]) / len(pair_spans)
        print(
            f"lmax = {lmax}, lmin = {lmin}, reads_included = {reads_included:.1%}",
        )
    return lmin, lmax


@dataclass
class Configs:
    MIN_MAPQ: int = 40
    MIN_BC_OVERLAP: int = 3
    BAM_FILE: str = None
    DEBUG: bool = False
    DIR: str = ""
    MAX_LINKED_DIST: int = 10_000
    NUM_THREADS: int = 1
    COMPRESSION_THREADS: int = 1
    SD_MULT: int = 2
    LMIN: int = None
    LMAX: int = None
    MIN_SV: int = None
    MIN_LEN: int = None
    MAX_LEN: int = 200000
    MIN_READS: int = 2
    MIN_DISCS: int = 2
    CANDIDATES: str = None
    BLACKLIST_FILE: str = None
    BLACKLIST: Dict[str, List[Tuple[int, int]]] = None

    def update(self, constants):
        if "min_mapq" in constants:
            self.MIN_MAPQ = int(constants["min_mapq"])

        if "k" in constants:
            self.MIN_BC_OVERLAP = int(constants["k"])

        if "bam_file" in constants and os.path.exists(constants["bam_file"]):
            self.BAM_FILE = constants["bam_file"]
        else:
            sys.exit(
                "Missing path to BAM file in config file. Please include this with the "
                "pattern 'bam_file=/path/to/file.bam'"
            )

        if "DEBUG" in constants:
            self.DEBUG = bool(constants["DEBUG"])

        if "outdir" in constants:
            self.DIR = constants["outdir"]
            if not os.path.exists(self.DIR):
                os.makedirs(self.DIR)

        if "d" in constants:
            self.MAX_LINKED_DIST = int(constants["d"])

        if "threads" in constants:
            self.NUM_THREADS = int(constants["threads"])

        if "sd_mult" in constants:
            self.SD_MULT = int(constants["sd_mult"])

        self.LMIN, self.LMAX = estimate_lmin_lmax(self.BAM_FILE, self.SD_MULT, self.DEBUG)

        if "min_sv" in constants:
            self.MIN_SV = int(constants["min_sv"])
        else:
            self.MIN_SV = 2 * self.LMAX

        if "min_len" in constants:
            self.MIN_LEN = int(constants["min_len"])
        else:
            self.MIN_LEN = 2 * self.LMAX

        if "max_len" in constants:
            self.MAX_LEN = int(constants["max_len"])

        if "min_reads" in constants:
            self.MIN_READS = int(constants["min_reads"])

        if "min_discs" in constants:
            self.MIN_DISCS = int(constants["min_discs"])

        if "candidates" in constants and os.path.exists(constants["candidates"]):
            self.CANDIDATES = constants["candidates"]

        if "blacklist" in constants and os.path.exists(constants["blacklist"]):
            self.BLACKLIST_FILE = constants["blacklist"]
            self.BLACKLIST = parse_blacklist(constants["blacklist"])

    @classmethod
    def from_file(cls, open_file):
        constants = cls.parse_configs(open_file)
        c = cls()
        c.update(constants)
        return c

    @staticmethod
    def parse_configs(iterator):
        fileconfigs = {}
        for line in iterator:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue

            name, val = line.split("=")
            name = name.strip(" ")
            val = val.strip(" ")
            fileconfigs[name] = val
        return fileconfigs
