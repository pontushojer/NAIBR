#!/usr/bin/env python3
"""
Convert NAIBR BEDPE files to VCF
"""
import argparse
import sys
import contextlib

from naibr.utils import NovelAdjacency, build_vcf_header


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


def get_chrom_lengths(reference):
    chrom_lengths = []
    with open(reference) as f:
        for line in f:
            chrom, length, *_ = line.strip().split("\t")
            assert str.isnumeric(length)
            chrom_lengths.append((chrom, int(length)))

    assert len(chrom_lengths) > 0, f"No chromosomes in reference '{reference}'"
    return chrom_lengths


def main():
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
    run(parser.parse_args())
    sys.exit(0)


def run(args):
    chrom_lengths = get_chrom_lengths(args.ref)

    vcf_header = build_vcf_header(chrom_lengths, sample=args.sample_name)

    with open(args.bedpe, "r") as bedpe_reader, smart_open(args.output_vcf) as vcf_writer:
        print(vcf_header, file=vcf_writer)
        nr = 0
        for entry in bedpe_reader:
            if entry.startswith("Chr1"):
                continue

            nr += 1

            novel_adjacency = NovelAdjacency.from_naibr(entry, add_chr=args.add_chr)
            vcf_entry = novel_adjacency.to_vcf(name=f"NAIBR_{nr:05}")
            print(vcf_entry, file=vcf_writer)


if __name__ == "__main__":
    main()
