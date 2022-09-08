from math import log
import os
from collections import defaultdict
import multiprocessing as mp
import itertools
import logging

import pysam

logger = logging.getLogger(__name__)


class NovelAdjacency:
    __slots__ = [
        "score",
        "score_by_hap",
        "chrm1",
        "break1",
        "chrm2",
        "break2",
        "orient",
        "haps",
        "discs",
        "pairs",
        "spans",
        "discs_by_hap",
        "pairs_by_hap",
        "spans_by_hap",
        "pass_threshold",
    ]

    def __init__(
        self,
        chrm1,
        chrm2,
        indi,
        indj,
        orient,
        score: float = None,
        haps=None,
        spans: int = None,
        discs: int = None,
        pairs: int = None,
        pass_threshold: bool = False,
    ):
        self.chrm1 = chrm1
        self.break1 = indi
        self.chrm2 = chrm2
        self.break2 = indj
        self.orient = orient

        self.score = -float("inf") if score is None else score
        self.score_by_hap = defaultdict(int)

        self.haps = (0, 0) if haps is None else haps
        self.spans = 0 if spans is None else spans
        self.discs = 0 if discs is None else discs
        self.pairs = 0 if pairs is None else pairs
        self.spans_by_hap = defaultdict(int)
        self.discs_by_hap = defaultdict(int)
        self.pairs_by_hap = defaultdict(int)

        self.pass_threshold = pass_threshold

    def __repr__(self):
        return (
            f"NovelAdjacency(chrm1={self.chrm1}, break1={self.break1}, chrm2={self.chrm2}, "
            f"break2={self.break2}, orient={self.orient}, haps={self.haps}, spans{self.spans}, "
            f"discs={self.discs}, pairs={self.pairs}, score={self.score}, "
            f"spans_by_hap={dict(self.spans_by_hap)}, discs_by_hap={dict(self.discs_by_hap)}, "
            f"pairs_by_hap={dict(self.pairs_by_hap)}, score_by_hap={dict(self.score_by_hap)}, "
            f"pass_threshold={self.pass_threshold}"
        )

    def __eq__(self, other):
        return self.coordinates() == other.coordinates()

    def __hash__(self):
        return hash(self.coordinates())

    def length(self):
        if not self.is_interchromosomal():
            return abs(self.break2 - self.break1)
        return 0

    def __len__(self):
        return self.length()

    def distance(self, other):
        if self.chrm1 == other.chrm1 and self.chrm2 == other.chrm2:
            return abs(self.break1 - other.break1), abs(self.break2 - other.break2)
        else:
            return float("inf"), float("inf")

    def coordinates(self):
        return self.chrm1, self.break1, self.chrm2, self.break2, self.orient

    def is_interchromosomal(self):
        return self.chrm1 != self.chrm2

    def add_spans(self, spans):
        for hap in spans:
            self.spans_by_hap[hap] += 1

    def add_score(self, candsplitmol):
        haplotype = candsplitmol.haplotype()
        self.score_by_hap[haplotype] += candsplitmol.score()
        self.pairs_by_hap[haplotype] += 1
        self.discs_by_hap[haplotype] += candsplitmol.discs

    def evaluate(self, thershold):
        if self.score > -float("inf"):
            self.pass_threshold = self.score >= thershold

    def svtype(self):
        """Translation of orientation to DEL, INV and DUP.
        Based on this: https://github.com/raphael-group/NAIBR/issues/11
        Note that this is not entirely accurate as the more complex variants are possible."""
        if self.is_interchromosomal():
            return "TRA"
        else:
            return {"+-": "DEL", "++": "INV", "--": "INV", "-+": "DUP"}.get(self.orient, "UNK")

    def homozygous(self) -> bool:
        """Return True if the variant is homozygous"""
        if self.haps == (3, 3) or 0 in set(self.haps):
            return True
        return False

    def genotype(self) -> str:
        """Return the genotype for the variant"""
        if self.homozygous():
            return "1/1"
        else:
            # TODO - Here we assume that if its the same haplotype on both sides of the
            #  NA it is phased. We need to know if the belong to the same phase set to
            #  acctually know this
            if self.haps == (1, 1):
                return "1|0"
            elif self.haps == (2, 2):
                return "0|1"
            else:  # haps either (1,2) or (2,1)
                return "0/1"

    def get_score(self):
        best_score = -float("inf")
        best_haps = (0, 0)
        total = 0
        for hap in list(self.score_by_hap):
            hap_score = self.score_by_hap[hap]

            if self.pairs_by_hap[hap] > self.spans_by_hap[hap]:
                total += hap_score

            if hap[0] != 0 and self.pairs_by_hap[(0, hap[1])] >= self.spans_by_hap[(0, hap[1])]:
                hap_score += max(0, self.score_by_hap[(0, hap[1])])

            if hap[1] != 0 and self.pairs_by_hap[(hap[0]), 0] >= self.spans_by_hap[(hap[0]), 0]:
                hap_score += max(0, self.score_by_hap[(hap[0]), 0])

            if hap[0] != 0 and hap[1] != 0 and self.pairs_by_hap[(0, 0)] >= self.spans_by_hap[(0, 0)]:
                hap_score += max(0, self.score_by_hap[(0, 0)])

            if (
                hap_score > best_score
                and self.pairs_by_hap[hap] > 0
                and self.pairs_by_hap[hap] >= self.spans_by_hap[hap]
            ):
                best_score = hap_score
                best_haps = hap

        pairs = sum(self.pairs_by_hap.values())
        spans = sum(self.spans_by_hap.values())
        discs = sum(self.discs_by_hap.values())

        # In this case there is good support from both haplotype pairs, (1,1) and (2,2)
        # or (1,2) and (2,1), indicating a homozygous call.
        if (
            -float("inf") < best_score < total  # TODO - why exclude cases where the best_score == total?
            and pairs >= spans
            and (
                (self.score_by_hap[(1, 1)] > 0 and self.score_by_hap[(2, 2)] > 0)
                or (self.score_by_hap[(1, 2)] > 0 and self.score_by_hap[(2, 1)] > 0)
            )
        ):
            self.haps = (3, 3)
            self.score = round(total, 3)
            self.discs = discs
            self.pairs = pairs
            self.spans = spans

        # In this case there is good support from one haplotype pair, indicating a heterozygous
        # call. Additional support from non-phased, haplotype=0, molecules counted towards call
        elif best_score > -float("inf"):
            self.haps = best_haps
            score = self.score_by_hap[best_haps]
            discs = self.discs_by_hap[best_haps]
            pairs = self.pairs_by_hap[best_haps]
            spans = self.spans_by_hap[best_haps]

            if best_haps[0] != 0 and self.score_by_hap[(0, best_haps[1])] > 0:
                score += self.score_by_hap[(0, best_haps[1])]
                discs += self.discs_by_hap[(0, best_haps[1])]
                pairs += self.pairs_by_hap[(0, best_haps[1])]
                spans += self.spans_by_hap[(0, best_haps[1])]

            if best_haps[1] != 0 and self.score_by_hap[(best_haps[0], 0)] > 0:
                score += self.score_by_hap[(best_haps[0], 0)]
                discs += self.discs_by_hap[(best_haps[0], 0)]
                pairs += self.pairs_by_hap[(best_haps[0], 0)]
                spans += self.spans_by_hap[(best_haps[0], 0)]

            if best_haps[1] != 0 and best_haps[0] != 0 and self.score_by_hap[(0, 0)] > 0:
                score += self.score_by_hap[(0, 0)]
                discs += self.discs_by_hap[(0, 0)]
                pairs += self.pairs_by_hap[(0, 0)]
                spans += self.spans_by_hap[(0, 0)]

            self.score = round(score, 3)
            self.discs = discs
            self.pairs = pairs
            self.spans = spans
        else:
            self.score = round(best_score, 3)
            self.discs = 0
            self.pairs = 0
            self.spans = 0

    @classmethod
    def from_naibr(cls, line: str, add_chr: bool = False):
        """Create a NovelAdjacency instance from a NAIBR-format BEDPE line. If add_chr
        is true 'chr' will be prepended to all chromosomes."""
        els = line.strip().split("\t")
        print(els[7])
        assert len(els) >= 10, "Less columns than expencted in BEDPE"
        return cls(
            chrm1=els[0] if not add_chr else f"chr{els[0]}",
            indi=int(els[1]),
            chrm2=els[2] if not add_chr else f"chr{els[2]}",
            indj=int(els[3]),
            orient=els[6],
            pairs=int(els[4]),
            discs=int(els[5]),
            haps=tuple(int(hap) for hap in els[7].split(",")),
            score=float(els[8]),
            pass_threshold=els[9] == "PASS",
        )

    def to_tuple(self):
        return (
            self.chrm1,
            self.break1,
            self.chrm2,
            self.break2,
            self.pairs,
            self.discs,
            self.orient,
            f"{self.haps[0]},{self.haps[1]}",
            self.score,
        )

    def to_naibr(self):
        """NAIBR default format"""
        return "\t".join(
            [
                self.chrm1,
                str(self.break1),
                self.chrm2,
                str(self.break2),
                str(self.pairs),
                str(self.discs),
                self.orient,
                f"{self.haps[0]},{self.haps[1]}",
                str(self.score),
                "PASS" if self.pass_threshold else "FAIL",
            ]
        )

    def to_bedpe(self, name: str):
        """10x-style BEDPE format according suitable for IGV visualization
        see https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bedpe"""
        info = ";".join(
            [
                f"LENGTH={len(self)}",
                f"NUM_SPLIT_MOLECULES={self.pairs}",
                f"NUM_DISCORDANT_READS={self.discs}",
                f"ORIENTATION={self.orient}",
                f"HAPLOTYPE={self.haps[0]},{self.haps[1]}",
                f"SVTYPE={self.svtype()}",
            ]
        )
        return "\t".join(
            [
                self.chrm1,
                str(self.break1),
                str(self.break1 + 1),
                self.chrm2,
                str(self.break2),
                str(self.break2 + 1),
                name,
                str(self.score),
                self.orient[0],  # Strand1
                self.orient[1],  # Strand 2
                "PASS" if self.pass_threshold else "FAIL",
                info,
            ]
        )

    def to_vcf(self, name):
        info = ";".join(
            (
                f"END={self.break2}",
                f"SVTYPE={self.svtype()}",
                f"SVLEN={len(self)}",
                f"SVSCORE={self.score}",
                f"NUM_SPLIT_MOLECULES={self.pairs}",
                f"NUM_DISCORDANT_READS={self.discs}",
                f"HAPLOTYPE={self.haps[0]},{self.haps[1]}",
            )
        )
        if self.is_interchromosomal():
            info += f";CHR1={self.chrm2}"

        return "\t".join(
            [
                self.chrm1,
                str(self.break1),
                name,
                "N",
                f"<{self.svtype()}>",
                str(self.score),
                "PASS" if self.pass_threshold else "FAIL",
                info,
                "GT",
                self.genotype(),
            ]
        )


class PERead:
    __slots__ = [
        "barcode",
        "orient",
        "chrm",
        "start",
        "end",
        "mean_mapq",
        "hap",
        "nextchrm",
        "nextstart",
        "nextend",
        "i",
        "j",
        "proper_pair",
    ]

    def __init__(self, read, barcode, mate=None):
        self.barcode = barcode
        self._set_read_orientation(read)

        self.chrm = read.reference_name
        self.start = read.reference_start
        self.end = read.reference_end
        self.hap = get_tag_default(read, "HP", default=0)

        self.nextchrm = read.next_reference_name
        self.nextstart = read.next_reference_start
        if mate is None:
            # Assume same aligment length and mapping quality for mate
            self.nextend = read.next_reference_start + read.reference_length
            self.mean_mapq = read.mapping_quality
        else:
            self.nextend = mate.reference_end
            self.mean_mapq = (mate.mapping_quality + read.mapping_quality) / 2

        self.i = self.start if read.is_reverse else self.end
        self.j = self.nextstart if read.mate_is_reverse else self.nextend
        self.proper_pair = read.is_proper_pair

    def __repr__(self):
        return (
            f"PERead(chrm={self.chrm}, start={self.start}, end={self.end}, nextchrm={self.nextchrm}, "
            f"nextstart={self.nextstart}, nextend={self.nextend}, barcode={self.barcode}, hap={self.hap}, "
            f"orient={self.orient}, i={self.i}, j={self.j}, proper_pair={self.proper_pair})"
        )

    def _set_read_orientation(self, read):
        a = "-" if read.is_reverse else "+"
        a += "-" if read.mate_is_reverse else "+"
        self.orient = a

    def is_discordant(self, min_sv):
        """Returns True if the read should be considered discordant"""
        # TODO - Which reads should acctually be considered discordant? Ususally it is base on orientation and
        #  pairwise distance.
        #
        # Note: The fact that MIN_SV is used as distance cutoff here is so we don't have candidate NAs
        # below this size.
        return (self.chrm != self.nextchrm or (self.j - self.i) > min_sv) and not self.proper_pair

    def is_concordant(self, lmin):
        # TODO - Which reads should acctually be considered concordant? From the paper it is stated: '
        #  A read-pair〈x,y〉 is concordant provided the distance between aligned reads f=ry−lx is
        #  between lmin and lmax and the orientations are ox=+,oy=-'. This is basically achived by
        #  `proper_pair` here. Why is to short fragments not concordant?
        return self.fragment_length() > lmin and self.proper_pair

    def fragment_length(self):
        return max(self.end, self.nextend) - min(self.start, self.nextstart)

    def mid(self):
        return int((self.i + self.j) / 2)


def first_read(read):
    chrom = read.reference_name
    mate_chrom = read.next_reference_name
    return (read.reference_start < read.next_reference_start and chrom == mate_chrom) or chrom < mate_chrom


def is_proper_chrom(chrom):
    return "Un" not in chrom and "random" not in chrom and "hs37d5" not in chrom


def get_read_orientation(read):
    a = "-" if read.is_reverse else "+"
    a += "-" if read.mate_is_reverse else "+"
    return a


def get_tag_default(read, tag, default=None):
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def plog(x, num):
    ret = 0
    for i in range(num):
        ret += safe_log(x)
    return ret


def safe_log(x):
    try:
        return log(x)
    except ValueError:
        return log(1e-300)


def is_close(pos1, pos2, max_dist):
    """True if read_pos is within lmax of break_pos"""
    return abs(pos1 - pos2) <= max_dist


def is_convergent(read, mate):
    return not read.is_reverse and mate.is_reverse


def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]


def threshold(cov):
    # In the publication the threshold relative coverage was determined by subsetting
    # down to 10X coverage. Below that it is uncertain how accurate it is.
    if cov < 10:
        logger.warning(f"Low coverage ({cov:.3f}X < 10X), the threshold value might not be accurate.")
    return round(6.943 * cov - 37.33, 3)


def filter_chrY(novel_adjacencies):
    filtered = []
    for na in novel_adjacencies:
        # TODO - Make this configurable to work with non-human genomes
        if "Y" not in na.chrm1 and "Y" not in na.chrm2:
            filtered.append(na)
    return filtered


def collapse_novel_adjacencies(novel_adjacencies, lmax):
    r = roundto(lmax, 100) * 5

    # Sort on decreasing score
    novel_adjacencies.sort(key=lambda x: x.score, reverse=True)

    # Keep NAs so that no other NA overlaps breakpoint positions within a distance of r
    collapsed_nas = {}
    for na in novel_adjacencies:
        already_appended = False
        norm_break1 = roundto(na.break1, r)
        norm_break2 = roundto(na.break2, r)

        # Check if already appended a novel_adjacency close to current
        for i, j in itertools.product([-r, 0, r], repeat=2):
            coords = (na.chrm1, norm_break1 + i, na.chrm2, norm_break2 + j)
            if coords in collapsed_nas:
                overlapping_na = collapsed_nas[coords]
                if abs(overlapping_na.break1 - na.break1) < r and abs(overlapping_na.break2 - na.break2) < r:
                    already_appended = True

        if not already_appended:
            collapsed_nas[(na.chrm1, norm_break1, na.chrm2, norm_break2)] = na

    return list(collapsed_nas.values())


def evaluate_threshold(novel_adjacencies, threshold_value):
    for na in novel_adjacencies:
        na.evaluate(threshold_value)


def get_chrom_lengths(bam_file):
    reads = pysam.AlignmentFile(bam_file, "rb")
    for chrom, length in zip(reads.references, reads.lengths):
        yield chrom, length


def build_vcf_header(chrom_lengths, sample=None):
    if sample is None:
        sample = "SAMPLE"
    """Build a VCF header. Requires a list of tuples with chromsome names and lengths"""
    header_string = """##fileformat=VCFv4.2
    ##source=NAIBR"""

    header_string += "".join([f"\n##contig=<ID={c},length={l}>" for c, l in chrom_lengths])

    header_string += f"""
    ##FILTER=<ID=PASS,Description="Passed the software filter">
    ##FILTER=<ID=FAIL,Description="Failed the software filter">
    ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, DUP=Duplication, INV=Inversion">
    ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
    ##INFO=<ID=SVSCORE,Number=1,Type=Float,Description="NAIBR Score for SV">
    ##INFO=<ID=NUM_DISCORDANT_READS,Number=1,Type=Integer,Description="Number of supporting discordant reads">
    ##INFO=<ID=NUM_SPLIT_MOLECULES,Number=1,Type=Integer,Description="Number of supporting split molecules">
    ##INFO=<ID=HAPLOTYPE,Number=1,Type=String,Description="Haplotype string from NAIBR">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype.">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{sample}"""
    return header_string


def write_novel_adjacencies(novel_adjacencies, directory, bam_file):
    fname = "NAIBR_SVs.bedpe"
    logger.info(f"Writing results to {os.path.join(directory, fname)}")
    with open(os.path.join(directory, fname), "w") as f:
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
        for na in novel_adjacencies:
            print(na.to_naibr(), file=f)

    fname2 = "NAIBR_SVs.reformat.bedpe"
    logger.info(f"Writing results to {os.path.join(directory, fname2)}")
    with open(os.path.join(directory, fname2), "w") as f:
        # This header is needed to recoginze it as a 10x-style BEDPE
        # see https://github.com/igvteam/igv/wiki/BedPE-Support
        print("#chrom1	start1	stop1	chrom2	start2	stop2	name	qual	strand1	strand2	filters	info", file=f)
        for nr, na in enumerate(novel_adjacencies, start=1):
            print(na.to_bedpe(f"NAIBR_{nr:05}"), file=f)

    fname3 = "NAIBR_SVs.vcf"
    logger.info(f"Writing results to {os.path.join(directory, fname3)}")
    with open(os.path.join(directory, fname3), "w") as f:
        header_string = build_vcf_header(get_chrom_lengths(bam_file))
        print(header_string, file=f)
        for nr, na in enumerate(novel_adjacencies, start=1):
            print(na.to_vcf(f"NAIBR_{nr:05}"), file=f)


def parallel_execute(function, input_list, threads=1):
    if threads > 1 and len(input_list) > 1:
        with mp.Pool(threads, maxtasksperchild=1) as pool:
            logger.info("running on %s threads" % str(threads))
            data = pool.map(function, input_list)
    else:
        data = map(function, input_list)
    return data


def roundto(number: int, step: int) -> int:
    """Round number to down the nearest multiple of step"""
    return number - number % step


class ConditionalFormatter(logging.Formatter):
    # Taken from https://stackoverflow.com/a/34954532
    def format(self, record):
        if hasattr(record, "nofmt") and record.nofmt:
            return record.getMessage()
        else:
            return logging.Formatter.format(self, record)


def input_candidates(openfile, min_sv: int = 1000):
    """Read input candidates from file. Filter out short candidates"""

    def parse_bedpe(lines):
        for line in lines:
            els = line.strip().split("\t")
            assert len(els) > 4
            assert els[1].isnumeric() and els[2].isnumeric() and els[-1] in {"+-",  "++", "--", "-+"}
            yield els[0], int(els[1]), els[3], int(els[4]), els[-1]

    def parse_naibr(lines):
        for line in lines:
            els = line.strip().split("\t")
            assert len(els) > 6
            assert els[1].isnumeric() and els[3].isnumeric() and els[6] in {"+-", "++", "--", "-+"}
            yield els[0], int(els[1]), els[2], int(els[3]), els[6]

    try:
        first_line = next(openfile)
    except StopIteration:
        first_line = ""

    cands = []
    if first_line.startswith("Chr1"):
        cands = list(parse_naibr(openfile))
    elif len(first_line.strip().split("\t")) > 4:
        openfile = itertools.chain([first_line], openfile)
        cands = list(parse_bedpe(openfile))

    n_total = len(cands)
    # Filter out short intrachromosomal candidates
    cands = [c for c in cands if abs(c[1] - c[3]) > min_sv or c[0] != c[2]]
    logger.info(f"Found {n_total:,} candidates in file of which {len(cands):,} are long enough")
    return cands
