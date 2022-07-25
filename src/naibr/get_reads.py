from __future__ import print_function, division
import pysam
import collections
import time
import copy

from .utils import *
from .global_vars import *
from .distributions import get_distributions


def inblacklist(cand):
    """
    Args: candidate novel adjacency
    Returns: True if candidate novel adjacency overlaps with
    interval in blacklist
    """
    if not BLACKLIST:
        return False
    for s, e in BLACKLIST[cand.chrm]:
        if s < cand.i < e:
            return True
    for s, e in BLACKLIST[cand.nextchrm]:
        if s < cand.j < e:
            return True
    return False


def make_barcodeDict(chrom):
    """
    Args: chromosome
    Function: Parses postion sorted BAM file and extracts relevant information
    Returns: dict barcoded reads, dict of reads overlapping ref positions,
    discordant reads, candidate NAs, total coverage
    """
    print("Getting candidates for chromosome %s" % chrom)
    cov = 0
    reads = pysam.AlignmentFile(BAM_FILE, "rb")
    reads_start = reads.lengths[list(reads.references).index(chrom)]
    reads_end = 0
    reads_by_LR = collections.defaultdict(list)
    discs_by_barcode = collections.defaultdict(list)
    discs = collections.defaultdict(list)
    interchrom_discs = collections.defaultdict(list)
    LRs_by_pos = collections.defaultdict(list)
    iterator = reads.fetch(chrom)
    t0 = time.time()
    nr = 0
    for nr, read in enumerate(iterator, start=1):
        reads_start = min(reads_start, read.reference_start)
        reads_end = max(reads_end, read.reference_end)
        cov += read.query_alignment_length
        ## DEBUG
        if DEBUG and cov > 100_000_000:
            break
        if pass_checks(read):
            peread = PERead(read)
            if peread.disc:
                peread.add_disc(read)
                discs_by_barcode[(peread.chrm, peread.nextchrm, peread.barcode)].append(peread)
                if peread.chrm == peread.nextchrm:
                    discs = add_disc(peread, discs)
                else:
                    interchrom_discs[
                        (
                            peread.chrm,
                            int(peread.i / LMAX) * LMAX,
                            peread.nextchrm,
                            int(peread.j / LMAX) * LMAX,
                            peread.orient,
                        )
                    ].append(peread)
            elif read.is_proper_pair and fragment_length(read) > LMIN:
                reads_by_LR[(peread.chrm, peread.barcode)].append(
                    (peread.start, peread.nextend, peread.hap, peread.mapq)
                )
                norm_mid = int(peread.mid() / MAX_LINKED_DIST) * MAX_LINKED_DIST
                if peread.barcode not in LRs_by_pos[(peread.chrm, norm_mid)]:
                    LRs_by_pos[(peread.chrm, norm_mid)].append(peread.barcode)

        if nr % 1_000_000 == 0:
            print(f"Processed {nr:,} reads ({time.time()-t0:.4f} s for last million).")
            t0 = time.time()

    # TODO - think about how to handle coverage a bit more accurate.
    #      - Use configuration to allow user input?
    #      - Otherwise calculate accurately across each chromosome (or subsection).
    coverage = cov / abs(reads_end-reads_start)
    print(f"Done reading chromosome {chrom}: coverage = {coverage:.3f}, reads = {nr:,}")
    return (
        reads_by_LR,
        LRs_by_pos,
        discs_by_barcode,
        discs,
        interchrom_discs,
        coverage,
    )


def signi(disc):
    if disc.orient[0] == "+":
        return disc.i + LMIN, disc.i + LMAX
    else:
        return disc.i - LMAX, disc.i - LMIN


def signj(disc):
    if disc.orient[1] == "+":
        return disc.j + LMIN, disc.j + LMAX
    else:
        return disc.j - LMAX, disc.j - LMIN


def add_disc(peread, discs):
    for a in [int(-LMAX / 2), 0, int(LMAX / 2)]:
        for b in [int(-LMAX / 2), 0, int(LMAX / 2)]:
            if (a == 0 or int((peread.i + a) / LMAX) * LMAX != int(peread.i / LMAX) * LMAX) and (
                b == 0 or int((peread.j + b) / LMAX) * LMAX != int(peread.j / LMAX) * LMAX
            ):
                discs[
                    (
                        peread.chrm,
                        int((peread.i + a) / LMAX) * LMAX,
                        peread.nextchrm,
                        int((peread.j + b) / LMAX) * LMAX,
                        peread.orient,
                    )
                ].append(peread)
    return discs


def largest_overlap(items):
    positions_i = collections.defaultdict(int)
    positions_j = collections.defaultdict(int)
    for disc in items:
        start, end = signi(disc)
        for pos in range(start, end):
            positions_i[pos] += 1
        start, end = signj(disc)
        for pos in range(start, end):
            positions_j[pos] += 1
    i = 0
    overlap_i = 0
    for pos, overlap in positions_i.items():
        if overlap > overlap_i:
            overlap_i = overlap
            i = pos
    j = 0
    overlap_j = 0
    for pos, overlap in positions_j.items():
        if overlap > overlap_j:
            overlap_j = overlap
            j = pos
    return i, j, overlap_i, overlap_j


def coordinates(s, e, orient):
    return s if orient == "+" else e


def disc_intersection(items):
    intersection = [0, float("Inf"), 0, float("Inf")]
    items.sort(key=lambda x: x.i)
    for disc in items:
        starti, endi = signi(disc)
        startj, endj = signj(disc)
        si, ei, sj, ej = intersection
        intersection = [max(si, starti), min(ei, endi), max(sj, startj), min(ej, endj)]
    si, ei, sj, ej = intersection
    if si and sj and si < ei and sj < ej:
        return si, ei, sj, ej
    return 0, 0, 0, 0


def get_candidates(discs, reads_by_LR):
    candidates = []
    p_len, p_rate, barcode_overlap = get_distributions(reads_by_LR)
    if p_len is None or p_rate is None:
        return None, None, None

    for key, items in discs.items():
        orient = key[4]
        si, ei, sj, ej = disc_intersection(items)
        if si and sj and len(items) >= MIN_DISCS:
            i = coordinates(si, ei, orient[0])
            j = coordinates(sj, ej, orient[1])
            cand = copy.copy(items[0])
            cand.i = i
            cand.j = j
            norm_i = int(cand.i / MAX_LINKED_DIST) * MAX_LINKED_DIST
            norm_j = int(cand.j / MAX_LINKED_DIST) * MAX_LINKED_DIST
            barcode_overlaps = barcode_overlap[(cand.chrm, norm_i, cand.nextchrm, norm_j)]
            if not inblacklist(cand) and (
                (cand.chrm == cand.nextchrm and cand.j - cand.i < MAX_LINKED_DIST)
                or barcode_overlaps >= MIN_BC_OVERLAP
            ):
                if not any(x.i == cand.i and x.j == cand.j for x in candidates):
                    candidates.append(cand)
    return candidates, p_len, p_rate


def make_barcodeDict_user(candidate):
    """
    Args: candidate novel adjacency input by user
    Function: Parses postion sorted BAM file and extracts positions indicated by candidate NAs
    Returns: dict barcoded reads, dict of reads overlapping ref positions,
    discordant reads, total coverage
    """
    w = MAX_LEN - MAX_LINKED_DIST
    cov = 0
    reads = pysam.AlignmentFile(BAM_FILE, "rb")
    reads_by_LR = collections.defaultdict(list)
    discs_by_barcode = collections.defaultdict(list)
    LRs_by_pos = collections.defaultdict(list)
    lengths = 0

    chrm1, break1, chrm2, break2, orientation = candidate
    if break1 > break2 or chrm1 > chrm2:
        chrm1, break1, chrm2, break2 = chrm2, break2, chrm1, break1
    if chrm1 != chrm2 or break1 + w < break2 - w:
        window = [(chrm1, break1 - w, break1 + w), (chrm2, break2 - w, break2 + w)]
    else:
        window = [(chrm1, break1 - w, break2 + w)]
    for chrm, s, e in window:
        lengths += e - s
        iterator = reads.fetch(chrm, max(0, s), e)
        for read in iterator:
            cov += read.query_alignment_length
            if pass_checks(read):
                peread = PERead(read)
                if peread.disc:
                    peread.add_disc(read)
                    discs_by_barcode[(peread.chrm, peread.nextchrm, peread.barcode)].append(peread)
                elif read.is_proper_pair and fragment_length(read) > LMIN:
                    reads_by_LR[(peread.chrm, peread.barcode)].append(
                        (peread.start, peread.nextend, peread.hap, peread.mapq)
                    )
                    norm_mid = int(peread.mid() / MAX_LINKED_DIST) * MAX_LINKED_DIST
                    if peread.barcode not in LRs_by_pos[(peread.chrm, norm_mid)]:
                        LRs_by_pos[(peread.chrm, norm_mid)].append(peread.barcode)

    cand = copy.copy(peread)
    cand.chrm = chrm1.strip("chr")
    cand.i = break1
    cand.nextchrm = chrm2.strip("chr")
    cand.j = break2
    cand.orient = orientation
    if (cand.chrm == cand.nextchrm and cand.i > cand.j) or (
        cand.chrm != cand.nextchrm and cand.chrm > cand.nextchrm
    ):
        cand.chrm = chrm2.strip("chr")
        cand.i = break2
        cand.nextchrm = chrm1.strip("chr")
        cand.j = break1
        cand.orient = orientation[::-1]

    coverage = cov / lengths
    if DEBUG:
        print(f"Candidate {candidate}: coverage = {coverage}")
    return reads_by_LR, LRs_by_pos, discs_by_barcode, [cand], coverage


def pass_checks(read):
    if read.mapping_quality < MIN_MAPQ:
        return False

    if read.is_duplicate:
        return False

    if read.is_secondary or read.is_supplementary:
        return False

    if not (is_proper_chrom(read.reference_name) and is_proper_chrom(read.next_reference_name)):
        return False

    if not read.has_tag("BX"):
        return False

    if not first_read(read):
        return False

    return True
