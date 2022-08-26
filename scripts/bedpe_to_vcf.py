"""
Convert NAIBR BEDPE files to VCF
"""
import argparse
import sys
import contextlib


@contextlib.contextmanager
def smart_open(filename=None):
    # From https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
    if filename and filename is not None and filename != "-":
        fh = open(filename, "w")
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


def create_header(sample_name, chrom_lengths):
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
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, All set to 1/1.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{sample_name}"""
    return header_string


def get_chrom_lengths(reference):
    chrom_lengths = []
    with open(reference) as f:
        for line in f:
            chrom, length, *_ = line.strip().split("\t")
            assert str.isnumeric(length)
            chrom_lengths.append((chrom, int(length)))

    assert len(chrom_lengths) > 0, f"No chromosomes in reference '{reference}'"
    return chrom_lengths


def get_svtype(orient, chrom1, chrom2):
    """Translation of orientation to DEL, INV and DUP.
    Based on this: https://github.com/raphael-group/NAIBR/issues/11
    Note that this is not entirely accurate as the more complex variants are possible."""
    if chrom1 != chrom2:
        return "TRA"
    else:
        return {"+-": "DEL", "++": "INV", "--": "INV", "-+": "DUP"}.get(orient, "UNK")


def main(args):
    chrom_lengths = get_chrom_lengths(args.ref)

    vcf_header = create_header(args.sample_name, chrom_lengths)

    with open(args.bedpe, "r") as bedpe_reader, smart_open(args.output_vcf) as vcf_writer:
        print(vcf_header, file=vcf_writer)
        nr = 0
        for entry in bedpe_reader:
            # Columns in BEDPE
            #    1. Chr1
            #    2. Break1
            #    3. Chr2
            #    4. Break2
            #    5. Split molecules
            #    6. Discordant reads
            #    7. Orientation
            #    8. Haplotype
            #    9. Score
            #    10. Pass filter

            if entry.startswith("Chr1"):
                continue

            nr += 1
            els = entry.strip().split("\t")
            assert len(els) >= 10, "Less columns than expencted in BEDPE"

            chrom1 = els[0] if not args.add_chr else f"chr{els[0]}"
            break1 = int(els[1])
            chrom2 = els[2] if not args.add_chr else f"chr{els[2]}"
            break2 = int(els[3])
            orient = els[6]
            pairs = float(els[4])
            discs = float(els[5])
            haps = els[7]
            score = float(els[8])
            pass_threshold = els[9] == "PASS"
            name = f"NAIBR_{nr:05}"
            length = break2 - break1 if chrom1 == chrom2 else 0
            svtype = get_svtype(orient, chrom1, chrom2)

            info = ";".join(
                (
                    f"END={break2}",
                    f"SVTYPE={svtype}",
                    f"SVLEN={length}",
                    f"SVSCORE={score}",
                    f"NUM_SPLIT_MOLECULES={pairs}",
                    f"NUM_DISCORDANT_READS={discs}",
                    f"HAPLOTYPE={haps}",
                )
            )
            if chrom1 != chrom2:
                info += f";CHR1={chrom2}"

            vcf_entry = "\t".join(
                [
                    chrom1,
                    str(break1),
                    name,
                    "N",
                    f"<{svtype}>",
                    str(score),
                    "PASS" if pass_threshold else "FAIL",
                    info,
                    "GT",
                    # TODO - This is not the correct genotype, but using './.' corresponds to
                    #  a "non-called" variant in downstream applications.
                    "1/1",
                ]
            )

            print(vcf_entry, file=vcf_writer)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    arg = parser.add_argument
    arg("bedpe", help="NAIBR-style BEDPE.")
    arg("-r", "--ref", required=True, help="List of chromosome lengths e.g. `*.fai`")
    arg("-o", "--output-vcf", default="-", help="Output VCF. Default: write to stdout")
    arg("--add-chr", action="store_true", help="Prepend 'chr' to chromsome names")
    arg(
        "-s",
        "--sample-name",
        type=lambda x: str(x).upper(),
        default="SAMPLE",
        help="Sample name. Default: %(default)s",
    )
    main(parser.parse_args())
