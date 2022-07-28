from __future__ import print_function, division
import math
import mpmath
import copy
import numpy as np
import linecache
import scipy.stats
import collections
import multiprocessing as mp

from .global_vars import *


class NovelAdjacency:
    __slots__ = [
        "score",
        "chrm1",
        "break1",
        "chrm2",
        "break2",
        "orient",
        "haps",
        "disc",
        "pairs",
        "spans",
    ]

    def __init__(self, chrm1, chrm2, indi, indj, orient):
        self.score = collections.defaultdict(int)
        self.chrm1 = chrm1
        self.break1 = indi
        self.chrm2 = chrm2
        self.break2 = indj
        self.orient = orient
        self.haps = (0, 0)
        self.spans = collections.defaultdict(int)
        self.disc = collections.defaultdict(int)
        self.pairs = collections.defaultdict(int)

    def add_spans(self, spans):
        for hap in spans:
            self.spans[hap] += 1

    def get_score(self):
        best_score = -float("inf")
        best_haps = (0, 0)
        total = 0
        for hap in list(self.score):
            s = self.score[hap]
            d = self.disc[hap]

            if self.pairs[hap] > self.spans[hap]:
                total += s

            if hap[0] != 0 and self.pairs[(0, hap[1])] >= self.spans[(0, hap[1])]:
                s += max(0, self.score[(0, hap[1])])
                d += self.disc[(0, hap[1])]

            if hap[1] != 0 and self.pairs[(hap[0]), 0] >= self.spans[(hap[0]), 0]:
                s += max(0, self.score[(hap[0]), 0])
                d += self.disc[(hap[0]), 0]

            if hap[0] != 0 and hap[1] != 0 and self.pairs[(0, 0)] >= self.spans[(0, 0)]:
                s += max(0, self.score[(0, 0)])
                d += self.disc[(0, 0)]

            if s > best_score and self.pairs[hap] > 0 and self.pairs[hap] >= self.spans[hap]:
                best_score = s
                best_haps = hap

        pairs = sum(self.pairs.values())
        spans = sum(self.spans.values())
        discs = sum(self.disc.values())

        if (
            -float("inf") < best_score < total
            and pairs >= spans
            and (
                (self.score[(1, 1)] > 0 and self.score[(2, 2)] > 0)
                or (self.score[(1, 2)] > 0 and self.score[(2, 1)] > 0)
            )
        ):
            self.haps = (3, 3)
            self.score = round(total, 3)
            self.disc = discs
            self.pairs = pairs
            self.spans = spans

        elif best_score > -float("inf"):
            self.haps = best_haps
            score = self.score[best_haps]
            discs = self.disc[best_haps]
            pairs = self.pairs[best_haps]
            spans = self.spans[best_haps]

            if best_haps[0] != 0 and self.score[(0, best_haps[1])] > 0:
                score += self.score[(0, best_haps[1])]
                discs += self.disc[(0, best_haps[1])]
                pairs += self.pairs[(0, best_haps[1])]
                spans += self.spans[(0, best_haps[1])]

            if best_haps[1] != 0 and self.score[(best_haps[0], 0)] > 0:
                score += self.score[(best_haps[0], 0)]
                discs += self.disc[(best_haps[0], 0)]
                pairs += self.pairs[(best_haps[0], 0)]
                spans += self.spans[(best_haps[0], 0)]

            if best_haps[1] != 0 and best_haps[0] != 0 and self.score[(0, 0)] > 0:
                score += self.score[(0, 0)]
                discs += self.disc[(0, 0)]
                pairs += self.pairs[(0, 0)]
                spans += self.spans[(0, 0)]

            self.score = round(score, 3)
            self.disc = discs
            self.pairs = pairs
            self.spans = spans
        else:
            self.score = round(best_score, 3)
            self.disc = 0
            self.pairs = 0
            self.spans = 0

    def to_tuple(self):
        return (
            self.chrm1,
            self.break1,
            self.chrm2,
            self.break2,
            self.pairs,
            self.disc,
            self.orient,
            f"{self.haps[0]},{self.haps[1]}",
            self.score,
        )


class LinkedRead:
    __slots__ = ["barcode", "chrm", "start", "disc", "end", "hap", "num"]

    def __init__(self, barcode):
        self.chrm = 0
        self.start = 0
        self.end = 0
        self.hap = 0
        self.num = 1
        self.disc = []
        self.barcode = barcode

    def add_num(self):
        self.num += 1

    def length(self):
        return float(self.end - self.start)


class PERead:
    __slots__ = [
        "barcode",
        "orient",
        "chrm",
        "start",
        "end",
        "mapq",
        "hap",
        "nextchrm",
        "nextstart",
        "nextend",
        "nextmapq",
        "nexthap",
        "i",
        "j",
        "disc",
    ]

    def __init__(self, read, barcode, mate=None):
        self.barcode = barcode
        self.orient = get_read_orientation(read)

        self.chrm = read.reference_name.strip("chr")
        self.start = read.reference_start
        self.end = read.reference_end
        self.mapq = read.mapping_quality
        self.hap = get_hap(read)

        self.nextchrm = read.next_reference_name.strip("chr")
        self.nextstart = read.next_reference_start
        if mate is None:
            # Assume same aligment length and mapping quality for mate
            self.nextend = read.next_reference_start + read.reference_length
            self.nextmapq = read.mapping_quality
        else:
            self.nextend = mate.reference_end
            self.nextmapq = mate.mapping_quality
        self.nexthap = self.hap

        self.i = self.start if read.is_reverse else self.end
        self.j = self.nextstart if read.mate_is_reverse else self.nextstart + read.reference_length
        self.j = self.nextstart if read.mate_is_reverse else self.nextend
        self.disc = self.is_disc() and not read.is_proper_pair

    def is_disc(self):
        return self.chrm != self.nextchrm or (self.j - self.i) > MIN_SV

    def fragment_length(self):
        return max(self.end, self.nextend) - min(self.start, self.nextstart)

    def mean_mapq(self):
        return int((self.mapq + self.nextmapq) / 2)

    def mid(self):
        if self.disc:
            return int((self.end + self.start) / 2)
        else:
            return int((self.i + self.j) / 2)


def first_read(read):
    chrom = read.reference_name
    mate_chrom = read.next_reference_name
    return (read.reference_start < read.next_reference_start and chrom == mate_chrom) or chrom < mate_chrom


def closest_linkedread(linkedread1, linkedread2, i):
    if not linkedread1:
        return linkedread2
    if not linkedread2:
        return linkedread1
    linkedread1_dist = min(abs(linkedread1[1] - i), abs(linkedread1[2]) - i)
    linkedread2_dist = min(abs(linkedread2[1] - i), abs(linkedread2[2]) - i)
    if linkedread1_dist < linkedread2_dist:
        return linkedread1
    return linkedread2


def is_proper_chrom(chrom):
    if chrom is not None:
        return "Un" not in chrom and "random" not in chrom and "hs37d5" not in chrom
    return False


def get_read_orientation(read):
    a = ""
    if read.is_reverse:
        a += "-"
    else:
        a += "+"
    if read.mate_is_reverse:
        a += "-"
    else:
        a += "+"
    return a


def get_hap(r):
    try:
        return r.get_tag("HP")
    except KeyError:
        return 0


def get_barcode(r):
    try:
        return r.get_tag("BX")
    except KeyError:
        return None


def plog(x, num):
    ret = 0
    for i in range(num):
        ret += log(x)
    return ret


def log(x):
    try:
        return math.log(x)
    except:
        return math.log(1e-300)


def is_close(index, read_pos, orient):
    if orient == "+":
        return abs(index - read_pos) <= LMAX
    else:
        return abs(index - read_pos) <= LMAX


def is_convergent(read, mate):
    return not read.is_reverse and mate.is_reverse


def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]


def threshold(cov):
    # In the publication the threshold relative coverage was determined by subsetting
    # down to 10X coverage. Below that it is uncertain how accurate it is.
    if cov < 10:
        print(f"WARNING: Low coverage ({cov:.3f}X < 10X), the threshold value might not be accurate.")
    return round(6.943 * cov - 37.33, 3)


def collapse(scores, threshold_value):
    def replace_23_24(chrom):
        replace = {"24": "Y", "23": "X"}
        return replace.get(chrom, chrom)

    l = []
    for linea in scores:
        chrom1, s1, chrom2, s2, split, disc, orient, haps, score = linea
        chrom1 = replace_23_24(chrom1)
        chrom2 = replace_23_24(chrom2)

        if "Y" not in chrom1 and "Y" not in chrom2:
            l.append(
                [chrom1, int(s1), chrom2, int(s2), float(split), float(disc), orient, haps, float(score)]
            )
    r = roundto(LMAX, 100) * 5

    # Sort on decreasing score
    l.sort(key=lambda x: x[-1], reverse=True)

    # Keep NAs so that no other NA overlaps breakpoint positions within a distance of r
    l2 = []
    nas = collections.defaultdict(list)
    for chrom1, s1, chrom2, s2, split, disc, orient, haps, score in l:
        already_appended = False
        norm_s1 = roundto(s1, r)
        norm_s2 = roundto(s2, r)
        for i in [-r, 0, r]:
            for j in [-r, 0, r]:
                if (chrom1, norm_s1 + i, chrom2, norm_s2 + j) in nas:
                    ch1, ds1, ch2, ds2 = nas[(chrom1, norm_s1 + i, chrom2, norm_s2 + j)][0:4]
                    if abs(ds1 - s1) < r and abs(ds2 - s2) < r:
                        already_appended = True
        if not already_appended:
            nas[(chrom1, norm_s1, chrom2, norm_s2)] = [
                chrom1,
                s1,
                chrom2,
                s2,
                split,
                disc,
                orient,
                haps,
                score,
            ]

    # Set PASS/FAIL based on threshold
    n_pass = 0
    for i, elem in nas.items():
        if elem[-1] >= threshold_value:
            l2.append(elem + ["PASS"])
            n_pass += 1
        else:
            l2.append(elem + ["FAIL"])

    if len(l2) > 0:
        print(f"Found {len(l2):,} SVs of which {n_pass:,} ({n_pass/len(l2):.1%}) passed the threshold")

    # Sort on descreasing score
    l2.sort(key=lambda x: (x[-2]), reverse=True)

    # Make list of lines for writing to file
    l2 = ["\t".join([str(y) for y in x]) + "\n" for x in l2]
    return l2


def write_scores(scores):
    fname = "NAIBR_SVs.bedpe"
    print(f"Writing results to {os.path.join(DIR, fname)}")
    with open(os.path.join(DIR, fname), "w") as f:
        print(
            "Chr1",
            "Break1",
            "Chr2",
            "Break2",
            "Split molecules",
            "Discordant reads",
            "Orientation",
            "Haplotype",
            "Score",
            "Pass filter",
            sep="\t",
            file=f,
        )
        for result in scores:
            f.write(result)


def parallel_execute(function, input_list):
    if NUM_THREADS != 1:
        pool = mp.Pool(NUM_THREADS, maxtasksperchild=1)
        map_fn = pool.map
        print("running on %s threads" % str(NUM_THREADS))
        data = map_fn(function, input_list)
        pool.close()
        pool.join()
    else:
        data = map(function, input_list)
    return data


def roundto(number: int, step: int) -> int:
    """Round number to down the nearset multiple of step"""
    return number - number % step
