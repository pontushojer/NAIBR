from __future__ import print_function, division
import sys
import time
import pysam
import collections

if len(sys.argv) < 2:
    sys.exit("No input config file")

from .get_reads import parse_candidate_region, get_distributions, parse_chromosome, get_candidates
from .global_vars import *
from .utils import flatten, parallel_execute, is_proper_chrom, write_scores
from .rank import predict_novel_adjacencies


def run_naibr_candidate(cand):
    """
    use user input candidate novel adjacencies
    """
    scores = 0
    reads_by_barcode, barcode_by_pos, discs_by_barcode, cands, coverage = parse_candidate_region(cand)
    p_len, p_rate, overlap = get_distributions(reads_by_barcode)
    if p_len is None:
        return scores
    scores = predict_novel_adjacencies(reads_by_barcode, barcode_by_pos, discs_by_barcode, cands, p_len, p_rate, coverage, False)
    return scores


def run_naibr(chrom):
    """
    automatically identify candidate novel adjacencies
    """
    (
        reads_by_barcode,
        barcode_by_pos,
        discs_by_barcode,
        discs,
        interchrom_discs,
        coverage,
    ) = parse_chromosome(chrom)
    scores = 0
    if len(reads_by_barcode) > 0:
        t_get_candidates = time.time()
        cands, p_len, p_rate = get_candidates(discs, reads_by_barcode)
        print(f"Got candidate noval adjacencies from data: {time.time() - t_get_candidates:.4f} s")
        if cands is None:
            print("No candidates from %s" % chrom)
            return (
                reads_by_barcode,
                barcode_by_pos,
                discs_by_barcode,
                interchrom_discs,
                coverage,
                scores,
            )

        print("ranking %i candidates from %s" % (len(cands), chrom))
        scores = predict_novel_adjacencies(
            reads_by_barcode,
            barcode_by_pos,
            discs_by_barcode,
            cands,
            p_len,
            p_rate,
            coverage,
            False,
        )
    else:
        print("No candidates from %s" % chrom)
    return reads_by_barcode, barcode_by_pos, discs_by_barcode, interchrom_discs, coverage, scores


def main():
    starttime = time.time()

    print("========= NAIBR =========")
    print("-------- CONFIGS --------")
    print("FILES")
    print(f"  bam_file = {BAM_FILE}")
    print(f"  candidates = {CANDIDATES}")
    print(f"  blacklist = {BLACKLIST_FILE}")
    print(f"  outdir = {DIR}")
    print("PARAMETERS")
    print(f"  d =         {MAX_LINKED_DIST:>9,} (maximum distance between read-pairs in a linked-read)")
    print(f"  min_mapq =  {MIN_MAPQ:>9}")
    print(f"  k =         {MIN_BC_OVERLAP:>9} (minimum number of barcode overlaps supporting a candidate NA)")
    print(f"  min_sv =    {MIN_SV:>9,} (minimum size of structural variant)")
    print(f"  sd_mult =   {SD_MULT:>9} (stddev multiplier for lmax/lmin calculcation)")
    print(f"  min_len =   {MIN_LEN:>9,}")
    print(f"  max_len =   {MAX_LEN:>9,}")
    print(f"  min_reads = {MIN_READS:>9}")
    print(f"  min_discs = {MIN_DISCS:>9}")
    print(f"  threads =   {NUM_THREADS:>9}")
    print(f"  DEBUG =     {str(DEBUG).rjust(9)}")
    print("-------------------------")

    if len(CANDIDATES) > 0:
        print("Using user defined candidates")
        cands = []
        with open(CANDIDATES) as f:
            for line in f:
                els = line.strip().split("\t")
                if len(els) > 4:
                    cands.append([els[0], int(els[1]), els[3], int(els[4]), els[-1]])
        scores = flatten(parallel_execute(run_naibr_candidate, cands))
        write_scores(scores)

    else:
        reads = pysam.AlignmentFile(BAM_FILE, "rb")
        chroms = reads.references
        chroms = [x for x in chroms if is_proper_chrom(x)]
        data = parallel_execute(run_naibr, chroms)
        reads_by_barcode = collections.defaultdict(list)
        barcode_by_pos = collections.defaultdict(list)
        discs_by_barcode = collections.defaultdict(list)
        discs = collections.defaultdict(list)
        coverage = []
        scores = []
        for (
            reads_by_barcode_chrom,
            barcode_by_pos_chrom,
            discs_by_barcode_chrom,
            discs_chrom,
            cov_chrom,
            scores_chrom,
        ) in data:
            if scores_chrom:
                reads_by_barcode.update(reads_by_barcode_chrom)
                barcode_by_pos.update(barcode_by_pos_chrom)
                discs_by_barcode.update(discs_by_barcode_chrom)
                discs.update(discs_chrom)
                coverage.append(cov_chrom)
                scores += scores_chrom
        cands, p_len, p_rate = get_candidates(discs, reads_by_barcode)
        if cands is not None:
            print("ranking %i interchromosomal candidates" % len(cands))
            scores += predict_novel_adjacencies(
                reads_by_barcode,
                barcode_by_pos,
                discs_by_barcode,
                cands,
                p_len,
                p_rate,
                np.mean(coverage),
                True,
            )
        else:
            print("No interchromosomal candidates")
        write_scores(scores)

    print("Finished in", (time.time() - starttime) / 60.0, "minutes")
    return


if __name__ == "__main__":
    main()
