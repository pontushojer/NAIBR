import collections
import logging
import os
from pathlib import Path
import sys
from dataclasses import dataclass, field
from typing import Mapping, Dict, List, Tuple, Optional, TYPE_CHECKING

import numpy as np
import pysam

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    PathType = os.PathLike[str]
else:
    # Remove when python >= 3.9 is required
    # see https://mypy.readthedocs.io/en/latest/runtime_troubles.html#generic-builtins
    PathType = os.PathLike


def parse_blacklist(fname: PathType):
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


def estimate_lmin_lmax(bam_file, sd_mult: int, debug: bool = False) -> Tuple[int, int]:
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

    pair_spans_array = np.array(pair_spans)
    mean_dist = pair_spans_array.mean()
    std_dist = pair_spans_array.std()

    lmin = max(int(mean_dist - std_dist * sd_mult), int(np.mean(reads_lengths)))
    lmax = max(int(mean_dist + std_dist * sd_mult), 100)

    if debug:
        # Estimate how large a percentage of reads have span within [lmin, lmax].
        reads_included = len(pair_spans_array[(pair_spans_array < lmax) & (pair_spans_array > lmin)]) / len(
            pair_spans_array
        )
        print(
            f"lmax = {lmax}, lmin = {lmin}, reads_included = {reads_included:.1%}",
        )
    return lmin, lmax


@dataclass
class Configs:
    MIN_MAPQ: int = 40
    MIN_BC_OVERLAP: int = 3
    BAM_FILE: Optional[os.PathLike[str]] = None
    DEBUG: bool = False
    DIR: os.PathLike = field(default=Path.cwd())
    MAX_LINKED_DIST: int = 10_000
    NUM_THREADS: int = 1
    COMPRESSION_THREADS: int = 1
    SD_MULT: int = 2
    LMIN: Optional[int] = None
    LMAX: Optional[int] = None
    MIN_SV: Optional[int] = None
    MIN_LEN: Optional[int] = None
    MAX_LEN: int = 200000
    MIN_READS: int = 2
    MIN_DISCS: int = 2
    CANDIDATES: Optional[os.PathLike[str]] = None
    BLACKLIST_FILE: Optional[os.PathLike[str]] = None
    BLACKLIST: Optional[Dict[str, List[Tuple[int, int]]]] = None
    _PREFIX: str = field(default="NAIBR_SVs")

    @property
    def PREFIX(self) -> str:
        return self._PREFIX

    @PREFIX.setter
    def PREFIX(self, value: str):
        if value == "":
            if self.BAM_FILE is not None:
                value = Path(self.BAM_FILE).stem
            else:
                value = self.PREFIX
        self._PREFIX = value

    def update(self, constants: Mapping[str, str]):
        if "min_mapq" in constants:
            self.MIN_MAPQ = int(constants["min_mapq"])

        if "k" in constants:
            self.MIN_BC_OVERLAP = int(constants["k"])

        if "bam_file" in constants and Path(constants["bam_file"]).exists():
            self.BAM_FILE = Path(constants["bam_file"])
        else:
            sys.exit(
                "Missing path to BAM file in config file. Please include this with the "
                "pattern 'bam_file=/path/to/file.bam'"
            )

        if "DEBUG" in constants:
            self.DEBUG = bool(constants["DEBUG"])

        if "outdir" in constants:
            self.DIR = Path(constants["outdir"])
            if self.DIR.exists() and not self.DIR.is_dir():
                sys.exit(f"{self.DIR} exists and is not a directory.")

            if not self.DIR.exists():
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

        if "candidates" in constants and Path(constants["candidates"]).exists():
            self.CANDIDATES = Path(constants["candidates"])

        if "blacklist" in constants and Path(constants["blacklist"]).exists():
            self.BLACKLIST_FILE = Path(constants["blacklist"])
            self.BLACKLIST = parse_blacklist(self.BLACKLIST_FILE)

        if "prefix" in constants:
            self.PREFIX = constants["prefix"]

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
