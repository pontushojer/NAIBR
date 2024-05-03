import collections
import logging
import sys
import time
from itertools import cycle

import numpy as np
import pysam

from .distributions import get_distributions, get_linkedread_distributions
from .utils import NovelAdjacency, PERead, get_tag_default, is_proper_chrom, roundto

logger = logging.getLogger(__name__)


def inblacklist(chrom, pos, blacklist):
    """
    True if position overlaps with interval in blacklist
    """
    if not blacklist:
        return False
    return any(s < pos < e for s, e in blacklist[chrom])


def parse_mapped_pairs(iterator, min_mapq):
    mates = {}
    for read in iterator:
        if read.mapping_quality < min_mapq:
            continue

        if read.is_duplicate or read.is_supplementary or read.is_secondary or not read.is_paired:
            continue

        if read.query_name in mates:
            yield mates.pop(read.query_name), read
        else:
            # Only check for unmapped mate as the any accessed reads here should be mapped.
            if read.mate_is_unmapped:
                continue

            if read.reference_name != read.next_reference_name:
                # Interchromosomal discordant reads
                if (
                    is_proper_chrom(read.next_reference_name)
                    and read.reference_name < read.next_reference_name
                ):
                    yield read, None

            mates[read.query_name] = read


def progress(iterator, desc=None, unit=None, interval=0.5):
    """Simple progress meter"""
    # Defaults
    desc = "Processed" if desc is None else desc
    unit = "els" if unit is None else unit

    # Initials
    t0 = time.perf_counter()
    t_start = time.perf_counter()
    nr = 0
    prev_nr = 0
    spinner = cycle(["-", "\\", "|", "/"])
    for nr, el in enumerate(iterator, start=1):
        yield el

        if time.perf_counter() - t0 > interval:
            step = nr - prev_nr
            prev_nr += step
            units_per_s = step // (time.perf_counter() - t0)
            print(
                f"\r{desc}: {next(spinner)} {nr:,} {unit} ({units_per_s:,.0f} {unit}/s)",
                file=sys.stderr,
                end="\r",
            )
            t0 = time.perf_counter()

    # Final printout
    print(
        f"\r{desc}: DONE! {nr:,} {unit} (total time = {time.perf_counter() - t_start:.3f} s)!",
        file=sys.stderr,
    )


def parse_chromosome(chrom, configs):
    """
    Args: chromosome
    Function: Parses postion sorted BAM file and extracts relevant information
    Returns: dict barcoded reads, dict of reads overlapping ref positions,
    discordant reads, candidate NAs, total coverage
    """
    logger.info("Getting candidates for chromosome %s" % chrom)
    cov = 0
    reads_by_barcode = collections.defaultdict(list)
    discs_by_barcode = collections.defaultdict(list)
    discs = collections.defaultdict(list)
    interchrom_discs = collections.defaultdict(list)
    barcodes_by_pos = collections.defaultdict(set)
    n_disc = 0
    n_conc = 0
    n_total = 0
    with pysam.AlignmentFile(configs.BAM_FILE, "rb", threads=configs.COMPRESSION_THREADS) as reads:
        iterator = reads.fetch(chrom)
        for read, mate in progress(
            parse_mapped_pairs(iterator, min_mapq=configs.MIN_MAPQ), desc=f"Processing {chrom}", unit="pairs"
        ):
            if mate is not None:
                cov += read.query_alignment_length + mate.query_alignment_length
            else:
                cov += read.query_alignment_length

            barcode = get_tag_default(read, "BX")
            if barcode is None:
                continue

            n_total += 1
            peread = PERead(read, barcode, mate=mate)
            if peread.is_discordant(min_sv=configs.MIN_SV):
                n_disc += 1
                discs_by_barcode[(peread.chrm, peread.nextchrm, peread.barcode)].append(peread)
                if peread.chrm == peread.nextchrm:
                    add_disc(peread, discs, lmax=configs.LMAX)
                else:
                    interchrom_discs[
                        (
                            peread.chrm,
                            roundto(peread.i, configs.LMAX),
                            peread.nextchrm,
                            roundto(peread.j, configs.LMAX),
                            peread.orient,
                        )
                    ].append(peread)
            elif peread.is_concordant(lmin=configs.LMIN):
                n_conc += 1
                reads_by_barcode[(peread.chrm, peread.barcode)].append(
                    (peread.start, peread.nextend, peread.hap, peread.mean_mapq)
                )
                norm_mid = roundto(peread.mid(), configs.MAX_LINKED_DIST)
                barcodes_by_pos[(peread.chrm, norm_mid)].add(peread.barcode)

    reads_start = float("inf")
    reads_end = 0
    readarray_by_barcode = {}
    dtype = [("start", int), ("end", int), ("hap", np.uint), ("mapq", int)]
    for barcode, reads in reads_by_barcode.items():
        readarray = np.array(reads, dtype=dtype)
        readarray.sort(order="start")
        readarray_by_barcode[barcode] = readarray
        reads_start = min(reads_start, readarray["start"][0])
        reads_end = max(reads_end, readarray["end"][-1])

    # TODO - think about how to handle coverage a bit more accurate.
    #      - Use configuration to allow user input?
    #      - Otherwise calculate accurately across each chromosome (or subsection).
    coverage = cov / (reads_end - reads_start)
    logger.info(f"Done reading chromosome {chrom}: coverage = {coverage:.3f}X")
    fraction_conc = n_conc / n_total if n_total != 0 else 0
    fraction_disc = n_disc / n_total if n_total != 0 else 0
    logger.info(
        f"{chrom}: total pairs: {n_total:,}, discordant: {n_disc:,} ({fraction_disc:.2%}), "
        f"concordant: {n_conc:,} ({fraction_conc:.2%})"
    )

    return (
        readarray_by_barcode,
        barcodes_by_pos,
        discs_by_barcode,
        discs,
        interchrom_discs,
        coverage,
    )


def signi(disc, configs):
    if disc.orient[0] == "+":
        return disc.i + configs.LMIN, disc.i + configs.LMAX
    else:
        return disc.i - configs.LMAX, disc.i - configs.LMIN


def signj(disc, configs):
    if disc.orient[1] == "+":
        return disc.j + configs.LMIN, disc.j + configs.LMAX
    else:
        return disc.j - configs.LMAX, disc.j - configs.LMIN


def add_disc(peread, discs, lmax):
    half_lmax = lmax // 2
    norm_i = roundto(peread.i, lmax)
    norm_j = roundto(peread.j, lmax)
    for a in [-half_lmax, 0, half_lmax]:
        norm_ia = roundto(peread.i + a, lmax)
        for b in [-half_lmax, 0, half_lmax]:
            norm_jb = roundto(peread.j + b, lmax)
            if (a == 0 or norm_ia != norm_i) and (b == 0 or norm_jb != norm_j):
                discs[(peread.chrm, norm_ia, peread.nextchrm, norm_jb, peread.orient)].append(peread)


def disc_intersection(disc_reads, configs):
    intersection = [0, float("inf"), 0, float("inf")]
    disc_reads.sort(key=lambda x: x.i)
    for disc_read in disc_reads:
        starti, endi = signi(disc_read, configs)
        startj, endj = signj(disc_read, configs)
        si, ei, sj, ej = intersection
        intersection = [max(si, starti), min(ei, endi), max(sj, startj), min(ej, endj)]
    si, ei, sj, ej = intersection
    if si and sj and si < ei and sj < ej:
        return si, ei, sj, ej
    return 0, 0, 0, 0


def get_candidates_from_discs(discs, barcode_overlap, configs):
    candidates = set()
    for position, disc_reads in progress(discs.items(), desc="Parsing discs", unit="pos"):
        # Skip positions with too few discs
        if len(disc_reads) < configs.MIN_DISCS:
            continue

        # Try to intersect read pair positions
        si, ei, sj, ej = disc_intersection(disc_reads, configs)
        if si == 0 or sj == 0:
            continue

        chrm = position[0]
        nextchrm = position[2]
        orient = position[4]
        i = si if orient[0] == "+" else ei
        j = sj if orient[1] == "+" else ej

        # Skip if positions overlaps blacklist
        if inblacklist(chrm, i, blacklist=configs.BLACKLIST) or inblacklist(
            nextchrm, j, blacklist=configs.BLACKLIST
        ):
            continue

        cand = NovelAdjacency(chrm1=chrm, chrm2=nextchrm, indi=i, indj=j, orient=orient)
        norm_i = roundto(i, configs.MAX_LINKED_DIST)
        norm_j = roundto(j, configs.MAX_LINKED_DIST)
        barcode_overlaps = barcode_overlap.get((chrm, norm_i, nextchrm, norm_j), 0)
        if (
            chrm == nextchrm and j - i < configs.MAX_LINKED_DIST
        ) or barcode_overlaps >= configs.MIN_BC_OVERLAP:
            candidates.add(cand)
    return list(candidates)


def get_candidates(discs, reads_by_barcode, configs):
    p_len, p_rate, barcode_overlap, barcode_linkedreads = get_distributions(reads_by_barcode, configs)
    if p_len is None or p_rate is None:
        return [], None, None, collections.defaultdict(list)

    candidates = get_candidates_from_discs(discs, barcode_overlap, configs)
    return candidates, p_len, p_rate, barcode_linkedreads


def get_interchrom_candidates(interchrom_discs, linkedreads_by_barcode, configs):
    p_len, p_rate, barcode_overlap = get_linkedread_distributions(linkedreads_by_barcode, configs)

    if p_len is None or p_rate is None:
        return [], None, None

    barcode_overlap = {
        position: barcodes for position, barcodes in barcode_overlap.items() if position[0] != position[2]
    }

    candidates = get_candidates_from_discs(interchrom_discs, barcode_overlap, configs)
    return candidates, p_len, p_rate


def parse_candidate_region(candidate, configs):
    """
    Args: candidate novel adjacency input by user
    Function: Parses postion sorted BAM file and extracts positions indicated by candidate NAs
    Returns: dict barcoded reads, dict of reads overlapping ref positions,
    discordant reads, total coverage
    """
    w = configs.MAX_LEN - configs.MAX_LINKED_DIST
    cov = 0
    reads_by_barcode = collections.defaultdict(list)
    discs_by_barcode = collections.defaultdict(list)
    barcodes_by_pos = collections.defaultdict(set)
    discs = collections.defaultdict(list)
    interchrom_discs = collections.defaultdict(list)
    lengths = 0

    chrm1, break1, chrm2, break2, orientation = candidate

    # Sort breakpoints if needed
    if (chrm1 == chrm2 and break1 > break2) or (chrm1 != chrm2 and chrm1 > chrm2):
        chrm1, break1, chrm2, break2 = chrm2, break2, chrm1, break1
        orientation = orientation[::-1]

    if chrm1 != chrm2 or break1 + w < break2 - w:
        window = [(chrm1, break1 - w, break1 + w), (chrm2, break2 - w, break2 + w)]
    else:
        window = [(chrm1, break1 - w, break2 + w)]

    with pysam.AlignmentFile(configs.BAM_FILE, "rb", threads=configs.COMPRESSION_THREADS) as reads:
        for chrm, s, e in window:
            lengths += e - s
            iterator = reads.fetch(chrm, max(0, s), e)
            for read, mate in progress(
                parse_mapped_pairs(iterator, min_mapq=configs.MIN_MAPQ),
                desc=f"Readning {chrm}:{s}-{e}",
                unit="pairs",
            ):
                if mate is not None:
                    cov += read.query_alignment_length + mate.query_alignment_length
                else:
                    cov += read.query_alignment_length

                barcode = get_tag_default(read, "BX")
                if barcode is None:
                    continue

                peread = PERead(read, barcode, mate=mate)
                if peread.is_discordant(min_sv=configs.MIN_SV):
                    discs_by_barcode[(peread.chrm, peread.nextchrm, peread.barcode)].append(peread)
                    if peread.chrm == peread.nextchrm:
                        add_disc(peread, discs, lmax=configs.LMAX)
                    else:
                        interchrom_discs[
                            (
                                peread.chrm,
                                roundto(peread.i, configs.LMAX),
                                peread.nextchrm,
                                roundto(peread.j, configs.LMAX),
                                peread.orient,
                            )
                        ].append(peread)
                elif peread.is_concordant(lmin=configs.LMIN):
                    reads_by_barcode[(peread.chrm, peread.barcode)].append(
                        (peread.start, peread.nextend, peread.hap, peread.mean_mapq)
                    )
                    norm_mid = roundto(peread.mid(), configs.MAX_LINKED_DIST)
                    barcodes_by_pos[(peread.chrm, norm_mid)].add(peread.barcode)

    readarray_by_barcode = {}
    dtype = [("start", int), ("end", int), ("hap", np.uint), ("mapq", int)]
    for barcode, reads in reads_by_barcode.items():
        readarray = np.array(reads, dtype=dtype)
        readarray.sort(order="start")
        readarray_by_barcode[barcode] = readarray

    cand = NovelAdjacency(chrm1=chrm1, chrm2=chrm2, indi=break1, indj=break2, orient=orientation)

    coverage = cov / lengths
    logger.debug(f"Candidate {candidate}: coverage = {coverage}")
    return readarray_by_barcode, barcodes_by_pos, discs_by_barcode, discs, cand, interchrom_discs, coverage
