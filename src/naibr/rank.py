from __future__ import print_function, division
import os
from functools import partial, cache
import math
import multiprocessing as mp
from .utils import log, NovelAdjacency, is_close, roundto, collapse, threshold
from .global_vars import *


class CandSplitMol:
    __slots__ = ["scoreHA", "scoreH0", "hap_i", "hap_j", "discs"]

    def __init__(self):
        self.scoreHA = 0
        self.scoreH0 = 0
        self.hap_i = 0
        self.hap_j = 0
        self.discs = 0

    def score(self, candidate_split, novel_adjacency, plen, prate):
        linkedread_i, discs_mapqs, linkedread_j = candidate_split
        # chrm,start,end,num,hap,barcode
        chrm_i, si, ei, mapq_i, hap_i, barcode_i = linkedread_i
        chrm_j, sj, ej, mapq_j, hap_j, barcode_j = linkedread_j
        num_i, num_j = len(mapq_i), len(mapq_j)
        hap_i = max(hap_i)
        hap_j = max(hap_j)
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
        scoreHA = log(
            plen(lenA) * prate(rateA) * mapqiA * mapqjA * mapqdiscA
            + plen(lenLi) * prate(rateLiA) * plen(lenLj) * prate(rateLjA) * mapqiA * mapqjA * mapqdiscA
        )
        if chrm_i == chrm_j:
            scoreH0 = log(
                plen(len0) * prate(rate0) * mapqi0 * mapqj0 * mapqdisc0
                + plen(lenLi) * prate(rateLi0) * plen(lenLj) * prate(rateLj0) * mapqi0 * mapqj0 * mapqdisc0
            )
        else:
            scoreH0 = log(
                0 + plen(lenLi) * prate(rateLi0) * plen(lenLj) * prate(rateLj0) * mapqi0 * mapqj0 * mapqdisc0
            )
        self.scoreHA = scoreHA
        self.scoreH0 = scoreH0
        self.hap_i = hap_i
        self.hap_j = hap_j
        self.discs = len(discs_mapqs)
        novel_adjacency.score[(hap_i, hap_j)] += self.scoreHA - self.scoreH0
        novel_adjacency.pairs[(hap_i, hap_j)] += 1
        novel_adjacency.disc[(hap_i, hap_j)] += self.discs

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


@cache
def phred_probability(mapq: int) -> int:
    return math.pow(10, (mapq / -10.0))


def score_pair(candidate, plen, prate, discs_by_barcode, barcodes_by_pos, reads_by_barcode):
    candidate_splits, spans = get_candidate_splits(
        candidate, barcodes_by_pos, reads_by_barcode, discs_by_barcode
    )
    novel_adjacency = NovelAdjacency(
        candidate.chrm, candidate.nextchrm, candidate.i, candidate.j, candidate.orient
    )
    for candidate_split in candidate_splits:
        c = CandSplitMol()
        c.score(candidate_split, novel_adjacency, plen, prate)
    for hap in spans:
        novel_adjacency.spans[hap] += 1
    return novel_adjacency, candidate


def near_i(x, candidate):
    chrm = x[0]
    start = x[1]
    end = x[2]
    neari = (
        (candidate.orient[0] == "+" and end - LMAX < candidate.i < end + MAX_LINKED_DIST)
        or (candidate.orient[0] == "-" and start + LMAX > candidate.i > start - MAX_LINKED_DIST)
    ) and candidate.chrm == chrm
    return neari


def near_j(x, candidate):
    chrm = x[0]
    start = x[1]
    end = x[2]
    nearj = (
        (candidate.orient[1] == "+" and end - LMAX < candidate.j < end + MAX_LINKED_DIST)
        or (candidate.orient[1] == "-" and start + LMAX > candidate.j > start - MAX_LINKED_DIST)
    ) and candidate.nextchrm == chrm
    return nearj


def spanning(x, candidate):
    start = x[1]
    end = x[2]
    span = start < candidate.i and end > candidate.j
    return span


def discs(candidate, barcode, discs_by_barcode):
    ds = discs_by_barcode[(candidate.chrm, candidate.nextchrm, barcode)]
    ds = [
        read.mean_mapq()
        for read in ds
        if candidate.chrm == read.chrm
        and candidate.nextchrm == read.nextchrm
        and is_close(candidate.i, read.i, read.orient[0])
        and is_close(candidate.j, read.j, read.orient[1])
        and candidate.orient == read.orient
    ]
    return ds


def crossing(start, end, i):
    return start < i < end


def linked_reads(barcode, chrm, candidate, reads_by_barcode):
    span = []
    reads = reads_by_barcode[(chrm, barcode)]
    reads.sort(key=lambda x: x[0])
    # chrm,start,end,num,hap,barcode
    current_linkedread = [0, 0, 0, [], [], 0]
    linkedreads = []
    for start, end, hap, mapq in reads:
        if current_linkedread[0] == 0 or start - current_linkedread[2] > MAX_LINKED_DIST:
            if (
                current_linkedread[0] != 0
                and len(current_linkedread[3]) >= MIN_READS
                and current_linkedread[2] - current_linkedread[1] >= MIN_LEN
            ):
                linkedreads.append(current_linkedread)
            current_linkedread = [chrm, start, end, [mapq], [hap], barcode]
        elif (
            candidate.chrm == candidate.nextchrm
            and current_linkedread[0] != 0
            and (current_linkedread[2] < candidate.i and start > candidate.j)
        ):
            linkedreads.append(current_linkedread)
            current_linkedread = [chrm, start, end, [mapq], [hap], barcode]
        else:
            current_linkedread[2] = max(current_linkedread[2], end)
            current_linkedread[3].append(mapq)
            current_linkedread[4].append(hap)
    if (
        current_linkedread[0] != 0
        and len(current_linkedread[3]) >= MIN_READS
        and current_linkedread[2] - current_linkedread[1] >= MIN_LEN
    ):
        linkedreads.append(current_linkedread)

    for lr1, lr2 in zip(linkedreads[:-1], linkedreads[1:]):
        if lr2[1] - lr1[2] < MAX_LINKED_DIST:
            lr1_hap = max(lr1[-2])
            lr2_hap = max(lr2[-2])
            span.append((lr1_hap, lr2_hap))
    return linkedreads, span


def get_linkedreads(candidate, barcodes, reads_by_barcode, discs_by_barcode):
    candidate_splits = []
    spans = []
    for barcode in barcodes:
        span = []
        linkedreads, s = linked_reads(barcode, candidate.chrm, candidate, reads_by_barcode)
        span += s
        if candidate.chrm != candidate.nextchrm:
            linkedreads2, s = linked_reads(barcode, candidate.nextchrm, candidate, reads_by_barcode)
            linkedreads += linkedreads2
        linkedreads_i, linkedreads_j = [[], []]
        for x in linkedreads:
            if near_i(x, candidate):
                linkedreads_i = x
            elif near_j(x, candidate):
                linkedreads_j = x
            if spanning(x, candidate):
                span += [(x[-2][0], x[-2][-1])]

        if linkedreads_i and linkedreads_j and linkedreads_i[1] < linkedreads_j[1]:
            discs_mapqs = discs(candidate, barcode, discs_by_barcode)
            candidate_split = [linkedreads_i, discs_mapqs, linkedreads_j]
            candidate_splits.append(candidate_split)
        spans += span
    return candidate_splits, spans


def get_candidate_splits(candidate, barcodes_by_pos, reads_by_barcode, discs_by_barcode):
    candidate_i_norm = roundto(candidate.i, MAX_LINKED_DIST)
    candidate_j_norm = roundto(candidate.j, MAX_LINKED_DIST)
    break1_barcodes = (
        barcodes_by_pos[(candidate.chrm, candidate_i_norm)]
        | barcodes_by_pos[(candidate.chrm, candidate_i_norm - MAX_LINKED_DIST)]
        | barcodes_by_pos[(candidate.chrm, candidate_i_norm + MAX_LINKED_DIST)]
    )
    break2_barcodes = (
        barcodes_by_pos[(candidate.nextchrm, candidate_j_norm)]
        | barcodes_by_pos[(candidate.nextchrm, candidate_j_norm - MAX_LINKED_DIST)]
        | barcodes_by_pos[(candidate.nextchrm, candidate_j_norm + MAX_LINKED_DIST)]
    )
    candidate_splits, spans = get_linkedreads(
        candidate, break1_barcodes.intersection(break2_barcodes), reads_by_barcode, discs_by_barcode
    )
    return candidate_splits, spans


def get_cand_score(
    candidates, is_interchrom, plen, prate, discs_by_barcode, barcodes_by_pos, reads_by_barcode
):
    scores = []
    score_pair_with_args = partial(
        score_pair,
        plen=plen,
        prate=prate,
        discs_by_barcode=discs_by_barcode,
        barcodes_by_pos=barcodes_by_pos,
        reads_by_barcode=reads_by_barcode,
    )
    if not is_interchrom or NUM_THREADS == 1:
        scores = map(score_pair_with_args, candidates)
    elif is_interchrom and NUM_THREADS != 1:
        pool = mp.Pool(min(5, NUM_THREADS), maxtasksperchild=1)
        scores = pool.map(score_pair_with_args, candidates)
        pool.close()
        pool.join()
    # TODO - Based on this if NUM_THREADS is 1 no interchrom candidates will be evaluated.

    rets = []
    for best, candidate in scores:
        best.get_score()
        if best.pairs > 0 and best.score > -float("inf"):
            ret = (
                best.chrm1,
                best.break1,
                best.chrm2,
                best.break2,
                best.pairs,
                best.disc,
                candidate.orient,
                f"{best.haps[0]},{best.haps[1]}",
                best.score,
            )
            rets.append(ret)
    return rets


def predict_novel_adjacencies(
    reads_by_barcode,
    barcodes_by_pos,
    discs_by_barcode,
    candidates,
    p_len,
    p_rate,
    cov,
    interchrom,
):
    plen = p_len
    prate = p_rate
    scores = get_cand_score(
        candidates, interchrom, plen, prate, discs_by_barcode, barcodes_by_pos, reads_by_barcode
    )
    scores = [x for x in scores if x]
    scores = collapse(scores, threshold(cov))
    return scores
