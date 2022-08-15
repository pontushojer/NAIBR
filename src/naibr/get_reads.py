import pysam
import collections
import time
import numpy as np

from .utils import PERead, roundto, is_proper_chrom, get_barcode, NovelAdjacency
from .global_vars import configs
from .distributions import get_distributions, get_linkedread_distributions


def inblacklist(chrom, pos):
    """
    True if position overlaps with interval in blacklist
    """
    if not configs.BLACKLIST:
        return False
    return any(s < pos < e for s, e in configs.BLACKLIST[chrom])


def parse_mapped_pairs(iterator):
    mates = {}
    for read in iterator:
        if read.mapping_quality < configs.MIN_MAPQ:
            continue

        if read.is_duplicate or read.is_supplementary or read.is_secondary:
            continue

        if read.query_name in mates:
            yield mates.pop(read.query_name), read
        else:
            if read.is_unmapped or read.mate_is_unmapped:
                continue

            if read.reference_name != read.next_reference_name:
                # Interchromosomal discordant reads
                if is_proper_chrom(read.next_reference_name) and read.reference_name < read.next_reference_name:
                    yield read, None

            mates[read.query_name] = read


def progress(iterator, desc=None, step=1_000_000, unit=None):
    """Simple progress meter"""
    desc = "Processed" if desc is None else desc
    unit = "els" if unit is None else unit
    nr = 0
    t0 = time.perf_counter()
    for nr, el in enumerate(iterator, start=1):
        yield el

        if nr % step == 0:
            print(f"{desc} {nr:,} {unit} ({1_000_000 // (time.perf_counter()-t0):,.0f} {unit}/s).")
            t0 = time.perf_counter()


def parse_chromosome(chrom):
    """
    Args: chromosome
    Function: Parses postion sorted BAM file and extracts relevant information
    Returns: dict barcoded reads, dict of reads overlapping ref positions,
    discordant reads, candidate NAs, total coverage
    """
    print("Getting candidates for chromosome %s" % chrom)
    cov = 0
    reads = pysam.AlignmentFile(configs.BAM_FILE, "rb", threads=configs.COMPRESSION_THREADS)

    reads_by_barcode = collections.defaultdict(list)
    discs_by_barcode = collections.defaultdict(list)
    discs = collections.defaultdict(list)
    interchrom_discs = collections.defaultdict(list)
    barcodes_by_pos = collections.defaultdict(set)
    iterator = reads.fetch(chrom)
    for read, mate in progress(parse_mapped_pairs(iterator), unit="pairs"):
        if mate is not None:
            cov += read.query_alignment_length + mate.query_alignment_length
        else:
            cov += read.query_alignment_length

        barcode = get_barcode(read)
        if barcode is None:
            continue

        peread = PERead(read, barcode, mate=mate)
        if peread.disc:
            discs_by_barcode[(peread.chrm, peread.nextchrm, peread.barcode)].append(peread)
            if peread.chrm == peread.nextchrm:
                add_disc(peread, discs)
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
        # TODO - consider removing minimum fragment_length as these are still properly
        #        paired with passing mapq.
        elif read.is_proper_pair and peread.fragment_length() > configs.LMIN:
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
    print(f"Done reading chromosome {chrom}: coverage = {coverage:.3f}")

    return (
        readarray_by_barcode,
        barcodes_by_pos,
        discs_by_barcode,
        discs,
        interchrom_discs,
        coverage,
    )


def signi(disc):
    if disc.orient[0] == "+":
        return disc.i + configs.LMIN, disc.i + configs.LMAX
    else:
        return disc.i - configs.LMAX, disc.i - configs.LMIN


def signj(disc):
    if disc.orient[1] == "+":
        return disc.j + configs.LMIN, disc.j + configs.LMAX
    else:
        return disc.j - configs.LMAX, disc.j - configs.LMIN


def add_disc(peread, discs):
    half_lmax = configs.LMAX / 2
    norm_i = roundto(peread.i, configs.LMAX)
    norm_j = roundto(peread.j, configs.LMAX)
    for a in [-half_lmax, 0, half_lmax]:
        norm_ia = roundto(peread.i + a, configs.LMAX)
        for b in [-half_lmax, 0, half_lmax]:
            norm_jb = roundto(peread.j + b, configs.LMAX)
            if (a == 0 or norm_ia != norm_i) and (b == 0 or norm_jb != norm_j):
                discs[(peread.chrm, norm_ia, peread.nextchrm, norm_jb, peread.orient)].append(peread)


def coordinates(s, e, orient):
    return s if orient == "+" else e


def disc_intersection(disc_reads):
    intersection = [0, float("Inf"), 0, float("Inf")]
    disc_reads.sort(key=lambda x: x.i)
    for disc_read in disc_reads:
        starti, endi = signi(disc_read)
        startj, endj = signj(disc_read)
        si, ei, sj, ej = intersection
        intersection = [max(si, starti), min(ei, endi), max(sj, startj), min(ej, endj)]
    si, ei, sj, ej = intersection
    if si and sj and si < ei and sj < ej:
        return si, ei, sj, ej
    return 0, 0, 0, 0


def get_candidates_from_discs(discs, barcode_overlap):
    candidates = []
    for position, disc_reads in discs.items():
        # Skip positions with too few discs
        if len(disc_reads) < configs.MIN_DISCS:
            continue

        # Try to intersect read pair positions
        si, ei, sj, ej = disc_intersection(disc_reads)
        if si == 0 or sj == 0:
            continue

        chrm = position[0]
        nextchrm = position[2]
        orient = position[4]
        i = coordinates(si, ei, orient[0])
        j = coordinates(sj, ej, orient[1])

        # Skip if positions overlaps blacklist
        if inblacklist(chrm, i) or inblacklist(nextchrm, j):
            continue

        cand = NovelAdjacency(chrm1=chrm, chrm2=nextchrm, indi=i, indj=j, orient=orient)
        norm_i = roundto(i, configs.MAX_LINKED_DIST)
        norm_j = roundto(j, configs.MAX_LINKED_DIST)
        barcode_overlaps = barcode_overlap[(chrm, norm_i, nextchrm, norm_j)]
        if (chrm == nextchrm and j - i < configs.MAX_LINKED_DIST) or barcode_overlaps >= configs.MIN_BC_OVERLAP:
            # Check that not already added
            if not any(cand == x for x in candidates):
                candidates.append(cand)
    return candidates


def get_candidates(discs, reads_by_barcode):
    p_len, p_rate, barcode_overlap, barcode_linkedreads = get_distributions(reads_by_barcode)
    if p_len is None or p_rate is None:
        return None, None, None, None

    candidates = get_candidates_from_discs(discs, barcode_overlap)
    return candidates, p_len, p_rate, barcode_linkedreads


def get_interchrom_candidates(interchrom_discs, linkedreads_by_barcode):
    p_len, p_rate, barcode_overlap = get_linkedread_distributions(linkedreads_by_barcode)

    if p_len is None or p_rate is None:
        return None, None, None

    candidates = get_candidates_from_discs(interchrom_discs, barcode_overlap)
    return candidates, p_len, p_rate


def parse_candidate_region(candidate):
    """
    Args: candidate novel adjacency input by user
    Function: Parses postion sorted BAM file and extracts positions indicated by candidate NAs
    Returns: dict barcoded reads, dict of reads overlapping ref positions,
    discordant reads, total coverage
    """
    w = configs.MAX_LEN - configs.MAX_LINKED_DIST
    cov = 0
    reads = pysam.AlignmentFile(configs.BAM_FILE, "rb", threads=configs.COMPRESSION_THREADS)
    reads_by_barcode = collections.defaultdict(list)
    discs_by_barcode = collections.defaultdict(list)
    barcodes_by_pos = collections.defaultdict(set)
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

    for chrm, s, e in window:
        lengths += e - s
        iterator = reads.fetch(chrm, max(0, s), e)
        for read, mate in progress(parse_mapped_pairs(iterator), unit="pairs"):
            if mate is not None:
                cov += read.query_alignment_length + mate.query_alignment_length
            else:
                cov += read.query_alignment_length

            barcode = get_barcode(read)
            if barcode is None:
                continue

            peread = PERead(read, barcode, mate=mate)
            if peread.disc:
                discs_by_barcode[(peread.chrm, peread.nextchrm, peread.barcode)].append(peread)
            # TODO - consider removing minimum fragment_length as these are still properly
            #        paired with passing mapq.
            elif read.is_proper_pair and peread.fragment_length() > configs.LMIN:
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
    if configs.DEBUG:
        print(f"Candidate {candidate}: coverage = {coverage}")
    return readarray_by_barcode, barcodes_by_pos, discs_by_barcode, [cand], coverage
