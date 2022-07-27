import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

BAM = Path("tests/data/example_chr21.bam").absolute()
CANDIDATES = Path("tests/data/example_candidates.bedpe").absolute()
BLACKLIST = Path("tests/data/example_blacklist.bedpe").absolute()


@dataclass
class Configs:
    """Write configs for NAIBR"""

    bam_file: os.PathLike
    outdir: os.PathLike
    min_mapq: int = 40
    blacklist: os.PathLike = None
    candidates: os.PathLike = None
    d: int = 10000
    min_sv: int = 1000
    threads: int = 1
    k: int = 3
    DEBUG: bool = True

    def write_to(self, file):
        with open(file, "w") as f:
            for k, v in vars(self).items():
                if v is None:
                    v = ""
                print(f"{k}={v}", file=f)


def same_pairwise_elements(list1, list2):
    assert len(list1) == len(list2)
    for l1, l2 in zip(list1, list2):
        assert l1 == l2


def test_candidates(tmp_path):
    configs = Configs(bam_file=BAM, candidates=CANDIDATES, outdir=tmp_path)
    configs.write_to(tmp_path / "naibr.configs")

    subprocess.run(["naibr", tmp_path / "naibr.configs"])
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "21	18759465	21	18906782	28.0	0.0	--	1,1	567.255	PASS".split("\t")
        )


def test_candidates_parallel(tmp_path):
    configs = Configs(bam_file=BAM, candidates=CANDIDATES, outdir=tmp_path, threads=2)
    configs.write_to(tmp_path / "naibr.configs")

    subprocess.run(["naibr", tmp_path / "naibr.configs"])
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "21	18759465	21	18906782	28.0	0.0	--	1,1	567.255	PASS".split("\t")
        )


def test_novel(tmp_path):
    configs = Configs(bam_file=BAM, outdir=tmp_path)
    configs.write_to(tmp_path / "naibr.configs")

    subprocess.run(["naibr", tmp_path / "naibr.configs"])
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "21	18842204	21	18860085	2.0	0.0	++	0,0	18.004	FAIL".split("\t")
        )


def test_novel_parallel(tmp_path):
    configs = Configs(bam_file=BAM, outdir=tmp_path, threads=2)
    configs.write_to(tmp_path / "naibr.configs")

    subprocess.run(["naibr", tmp_path / "naibr.configs"])
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "21	18842204	21	18860085	2.0	0.0	++	0,0	18.004	FAIL".split("\t")
        )


def test_consistent(tmp_path):
    # Run without candidates
    configs = Configs(bam_file=BAM, outdir=tmp_path / "novel")
    configs.write_to(tmp_path / "novel.configs")

    subprocess.run(["naibr", tmp_path / "novel.configs"])
    output = tmp_path / "novel" / "NAIBR_SVs.bedpe"

    # Generate candidates.bedpe base on results
    new_candidates = tmp_path / "candidates.bedpe"
    with open(new_candidates, "w") as outfile, open(output) as infile:
        line = infile.readlines()[1]
        chr1, pos1, chr2, pos2, _, _, orient, *_ = line.split("\t")
        print(f"chr{chr1}", pos1, pos1, f"chr{chr2}", pos2, pos2, orient, sep="\t", file=outfile)

    # Run with generated candidates
    configs = Configs(bam_file=BAM, candidates=new_candidates, outdir=tmp_path / "rerun")
    configs.write_to(tmp_path / "rerun.configs")

    subprocess.run(["naibr", tmp_path / "rerun.configs"])
    output_rerun = tmp_path / "rerun" / "NAIBR_SVs.bedpe"
    with open(output_rerun) as f:
        new_line = f.readlines()[1]
        same_pairwise_elements(line, new_line)
