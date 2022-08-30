import collections
from functools import partial, lru_cache
import math
import numpy as np
import logging

from .utils import (
    safe_log,
    NovelAdjacency,
    is_close,
    roundto,
    filter_chrY,
    collapse_novel_adjacencies,
    threshold,
    evaluate_threshold,
)

logger = logging.getLogger(__name__)


class CandSplitMol:
    __slots__ = ["scoreHA", "scoreH0", "hap_i", "hap_j", "discs"]

    def __init__(self):
        self.scoreHA = 0
        self.scoreH0 = 0
        self.hap_i = 0
        self.hap_j = 0
        self.discs = 0

    def calculate_score(self, candidate_split, plen, prate):
        linkedread_i, discs_mapqs, linkedread_j = candidate_split
        # chrm,start,end,num,hap,barcode
        chrm_i, si, ei, mapq_i, hap_i, barcode_i = linkedread_i
        chrm_j, sj, ej, mapq_j, hap_j, barcode_j = linkedread_j
        num_i, num_j = len(mapq_i), len(mapq_j)
        hap_i = majority_vote(hap_i)
        hap_j = majority_vote(hap_j)
        len0 = ej - si
        rate0 = (num_j + num_i) / float(len0)
        lenA = (ei - si) + (ej - sj)
        rateA = (num_j + num_i) / float(lenA)
        lenLi = ei - si
        rateLi0 = (num_i) / float(lenLi)
        rateLiA = (num_i) / float(lenLi)
        lenLj = ej - sj
        rateLj0 = (num_j) / float(lenLj)
        rateLjA = (num_j) / float(lenLj)
        mapqiA = self.map_prob(mapq_i)
        mapqi0 = self.map_prob(mapq_i)
        mapqjA = self.map_prob(mapq_j)
        mapqj0 = self.map_prob(mapq_j)
        mapqdiscA = self.map_prob(discs_mapqs)
        mapqdisc0 = self.mismap_prob(discs_mapqs)
        scoreHA = safe_log(
            plen(lenA) * prate(rateA) * mapqiA * mapqjA * mapqdiscA
            + plen(lenLi) * prate(rateLiA) * plen(lenLj) * prate(rateLjA) * mapqiA * mapqjA * mapqdiscA
        )
        if chrm_i == chrm_j:
            scoreH0 = safe_log(
                plen(len0) * prate(rate0) * mapqi0 * mapqj0 * mapqdisc0
                + plen(lenLi) * prate(rateLi0) * plen(lenLj) * prate(rateLj0) * mapqi0 * mapqj0 * mapqdisc0
            )
        else:
            scoreH0 = safe_log(
                0 + plen(lenLi) * prate(rateLi0) * plen(lenLj) * prate(rateLj0) * mapqi0 * mapqj0 * mapqdisc0
            )
        self.scoreHA = scoreHA
        self.scoreH0 = scoreH0
        self.hap_i = hap_i
        self.hap_j = hap_j
        self.discs = len(discs_mapqs)

    def score(self):
        return self.scoreHA - self.scoreH0

    def haplotype(self):
        return self.hap_i, self.hap_j

    def map_prob(self, mapqs):
        prob = 1
        for mapq in mapqs:
            prob *= 1 - phred_probability(mapq)
        return prob

    def mismap_prob(self, mapqs):
        prob = 1
        for mapq in mapqs:
            prob *= phred_probability(mapq)
        return prob


@lru_cache(maxsize=None)
def phred_probability(mapq: int) -> float:
    return math.pow(10, (mapq / -10.0))


def near_i(x, candidate, configs):
    chrm = x[0]
    start = x[1]
    end = x[2]
    neari = (
        (candidate.orient[0] == "+" and end - configs.LMAX < candidate.break1 < end + configs.MAX_LINKED_DIST)
        or (
            candidate.orient[0] == "-"
            and start + configs.LMAX > candidate.break1 > start - configs.MAX_LINKED_DIST
        )
    ) and candidate.chrm1 == chrm
    return neari


def near_j(x, candidate, configs):
    chrm = x[0]
    start = x[1]
    end = x[2]
    nearj = (
        (candidate.orient[1] == "+" and end - configs.LMAX < candidate.break2 < end + configs.MAX_LINKED_DIST)
        or (
            candidate.orient[1] == "-"
            and start + configs.LMAX > candidate.break2 > start - configs.MAX_LINKED_DIST
        )
    ) and candidate.chrm2 == chrm
    return nearj


def spanning(linkedread, candidate):
    """Return True if the candidate is within the bounds of the linked reads"""
    start = linkedread[1]
    end = linkedread[2]
    span = start < candidate.break1 and end > candidate.break2
    return span


def discs(candidate, barcode, discs_by_barcode, lmax):
    ds = discs_by_barcode[(candidate.chrm1, candidate.chrm2, barcode)]
    ds = [
        read.mean_mapq
        for read in ds
        if candidate.chrm1 == read.chrm
        and not candidate.is_interchromosomal()
        and is_close(candidate.break1, read.i, read.orient[0], lmax)
        and is_close(candidate.break2, read.j, read.orient[1], lmax)
        and candidate.orient == read.orient
    ]
    return ds


def majority_vote(haps):
    """From a list of haplotypes select the most numerous that's not 0 (unphased)
    unless all are 0"""
    counts = collections.Counter(haps)
    if counts[1] == 0 and counts[2] == 0:
        return 0
    elif counts[1] > counts[2]:
        return 1
    else:
        return 2


def linked_reads(barcode, chrm, candidate, reads_by_barcode, configs):
    if (chrm, barcode) not in reads_by_barcode:
        return [], []

    reads = reads_by_barcode[(chrm, barcode)]

    # Select only reads with proximal to the start and end of the candidate
    gt_d_start = reads["start"] > (candidate.break1 - 2 * configs.MAX_LINKED_DIST)
    lt_d_end = reads["end"] < (candidate.break2 + 2 * configs.MAX_LINKED_DIST)

    proximal_reads = reads[gt_d_start & lt_d_end].copy()
    if len(proximal_reads) == 0:
        return [], []

    # Calculate the distance between neighbouring reads
    pair_dists = proximal_reads["start"][1:] - proximal_reads["end"][:-1]

    # Get indexes where the distance of neigbouring reads exceeds the MAX_LINKED_DIST.
    # These reads most likely originate from different molecuels.
    dist_break = pair_dists > configs.MAX_LINKED_DIST

    # Get indexes neigbouring reads span the candidate. These reads could originate the
    # same molecule. Not for interchromsomal candidates.
    cand_break = np.full(np.shape(dist_break), False)
    if not candidate.is_interchromosomal():
        cand_break = (proximal_reads["end"][:-1] < candidate.break1) & (
            proximal_reads["start"][1:] > candidate.break2
        )

    breaks = np.where(dist_break | cand_break)[0] + 1

    # List whether breaks originate only from reads spanning the candidate
    # These are not filtered the same way as others.
    break_only_cands = [b not in set(np.where(dist_break)[0] + 1) for b in breaks] + [False]
    linkedreads = []

    # If the prevous reads were split because of spanning the candidate the next reads should
    # be filtered the same way
    prev_break_only_cands = False
    for reads_group, only_cands in zip(np.split(proximal_reads, breaks), break_only_cands):
        start = reads_group["start"].min()
        end = reads_group["end"].max()
        mapq = list(reads_group["mapq"])
        hap = list(reads_group["hap"])
        nr_reads = len(mapq)
        if (only_cands or prev_break_only_cands) or (
            nr_reads >= configs.MIN_READS and end - start >= configs.MIN_LEN
        ):
            linkedreads.append([chrm, start, end, mapq, hap, barcode])

        prev_break_only_cands = only_cands

    span = []
    if not candidate.is_interchromosomal():
        for lr1, lr2 in zip(linkedreads[:-1], linkedreads[1:]):
            if lr2[1] - lr1[2] < configs.MAX_LINKED_DIST:
                lr1_hap = majority_vote(lr1[-2])
                lr2_hap = majority_vote(lr2[-2])
                span.append((lr1_hap, lr2_hap))
    return linkedreads, span


def get_linkedreads(candidate, barcodes, reads_by_barcode, discs_by_barcode, is_interchrom, configs):
    candidate_splits = []
    spans = []
    for barcode in barcodes:
        span = []
        if is_interchrom:
            # When is_interchrom reads_by_barcode is a map of chromosomes to barcodes to linkedreads.
            linkedreads = reads_by_barcode[candidate.chrm1][barcode]
            linkedreads2 = reads_by_barcode[candidate.chrm2][barcode]
            linkedreads.extend(linkedreads2)
        else:
            linkedreads, s = linked_reads(barcode, candidate.chrm1, candidate, reads_by_barcode, configs)
            span.extend(s)

        linkedread_i, linkedread_j = [None, None]
        for linkedread in linkedreads:
            if near_i(linkedread, candidate, configs):
                linkedread_i = linkedread
            elif near_j(linkedread, candidate, configs):
                linkedread_j = linkedread

            # A single linkedread cannot be spanning two chromosomes
            if not is_interchrom and spanning(linkedread, candidate):
                span += [(linkedread[-2][0], linkedread[-2][-1])]  # First and last haplotype on linked read

        if linkedread_i and linkedread_j and linkedread_i[1] < linkedread_j[1]:
            discs_mapqs = discs(candidate, barcode, discs_by_barcode, lmax=configs.LMAX)
            candidate_split = [linkedread_i, discs_mapqs, linkedread_j]
            candidate_splits.append(candidate_split)
        spans.extend(span)
    candidate.add_spans(spans)
    return candidate_splits


def get_candidate_splits(
    candidate: NovelAdjacency,
    barcodes_by_pos,
    reads_by_barcode,
    discs_by_barcode,
    is_interchrom,
    configs,
):
    # Get normalized positions for candidate split molecule
    candidate_break1_norm = roundto(candidate.break1, configs.MAX_LINKED_DIST)
    candidate_break2_norm = roundto(candidate.break2, configs.MAX_LINKED_DIST)

    # Create sets of proximal barcodes
    break1_barcodes = (
        barcodes_by_pos[(candidate.chrm1, candidate_break1_norm)]
        | barcodes_by_pos[(candidate.chrm1, candidate_break1_norm - configs.MAX_LINKED_DIST)]
        | barcodes_by_pos[(candidate.chrm1, candidate_break1_norm + configs.MAX_LINKED_DIST)]
    )
    break2_barcodes = (
        barcodes_by_pos[(candidate.chrm2, candidate_break2_norm)]
        | barcodes_by_pos[(candidate.chrm2, candidate_break2_norm - configs.MAX_LINKED_DIST)]
        | barcodes_by_pos[(candidate.chrm2, candidate_break2_norm + configs.MAX_LINKED_DIST)]
    )

    barcodes_intersect = break1_barcodes.intersection(break2_barcodes)
    candidate_splits = get_linkedreads(
        candidate,
        barcodes_intersect,
        reads_by_barcode,
        discs_by_barcode,
        is_interchrom,
        configs,
    )
    return candidate_splits


def score_pair(
    candidate, plen, prate, discs_by_barcode, barcodes_by_pos, reads_by_barcode, is_interchrom, configs
):
    candidate_splits = get_candidate_splits(
        candidate, barcodes_by_pos, reads_by_barcode, discs_by_barcode, is_interchrom, configs
    )
    for candidate_split in candidate_splits:
        c = CandSplitMol()
        c.calculate_score(candidate_split, plen, prate)
        candidate.add_score(c)

    return candidate


def get_cand_score(
    candidates, is_interchrom, plen, prate, discs_by_barcode, barcodes_by_pos, reads_by_barcode, configs
):
    score_pair_with_args = partial(
        score_pair,
        plen=plen,
        prate=prate,
        discs_by_barcode=discs_by_barcode,
        barcodes_by_pos=barcodes_by_pos,
        reads_by_barcode=reads_by_barcode,
        is_interchrom=is_interchrom,
        configs=configs,
    )
    scores = map(score_pair_with_args, candidates)

    rets = []
    for novel_adjacency in scores:
        novel_adjacency.get_score()
        if novel_adjacency.pairs > 0 and novel_adjacency.score > -float("inf"):
            rets.append(novel_adjacency)
    return rets


def predict_novel_adjacencies(
    reads_by_barcode,
    barcodes_by_pos,
    discs_by_barcode,
    candidates,
    plen,
    prate,
    cov,
    interchrom,
    configs,
):
    novel_adjacencies = get_cand_score(
        candidates, interchrom, plen, prate, discs_by_barcode, barcodes_by_pos, reads_by_barcode, configs
    )
    novel_adjacencies = filter_chrY(novel_adjacencies)
    novel_adjacencies = collapse_novel_adjacencies(novel_adjacencies, lmax=configs.LMAX)
    evaluate_threshold(novel_adjacencies, threshold(cov))
    return novel_adjacencies
