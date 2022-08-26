import io
from pathlib import Path

import pytest

from naibr.__main__ import main, run
from naibr.global_vars import Configs

BAM = Path("tests/data/example_chr21.bam").absolute()
CANDIDATES = Path("tests/data/example_candidates.bedpe").absolute()
CANDIDATES_REV = Path("tests/data/example_candidates_reversed.bedpe").absolute()
BLACKLIST = Path("tests/data/example_blacklist.bed").absolute()


def same_pairwise_elements(list1, list2):
    assert len(list1) == len(list2)
    for l1, l2 in zip(list1, list2):
        assert l1 == l2


def test_help_return_zero():
    for arg in ["-h", "--help", "help"]:
        with pytest.raises(SystemExit):
            errorcode = main([arg])
            assert errorcode == 0


def test_empty_return_zero():
    with pytest.raises(SystemExit):
        errorcode = main([])
        assert errorcode == 0


def test_missing_config_raises_error():
    with pytest.raises(SystemExit):
        errorcode = main(["some_non_existing.config"])
        assert errorcode != 0


def test_missing_bam_file_raises_error(tmp_path):
    config_file_without_bam = io.StringIO("")
    with pytest.raises(SystemExit):
        _ = Configs.from_file(config_file_without_bam)


def test_candidates(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM}\ncandidates={CANDIDATES}\noutdir={tmp_path}\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "chr21	18759465	chr21	18906782	28	0	--	1,1	567.174	PASS".split("\t"),
        )


def test_blacklist(tmp_path):
    assert BLACKLIST.exists()
    config_file = io.StringIO(f"bam_file={BAM}\nblacklist={BLACKLIST}\noutdir={tmp_path}\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "chr21	18842204	chr21	18860085	2	0	++	0,0	18.002	FAIL".split("\t"),
        )


def test_candidates_reversed(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM}\ncandidates={CANDIDATES_REV}\noutdir={tmp_path}\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "chr21	18759465	chr21	18906782	28	0	--	1,1	567.174	PASS".split("\t"),
        )


def test_candidates_parallel(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM}\ncandidates={CANDIDATES}\noutdir={tmp_path}\nthreads=2\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "chr21	18759465	chr21	18906782	28	0	--	1,1	567.174	PASS".split("\t"),
        )


def test_novel(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM}\noutdir={tmp_path}\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "chr21	18842204	chr21	18860085	2	0	++	0,0	18.002	FAIL".split("\t"),
        )


def test_novel_parallel(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM}\noutdir={tmp_path}\nthreads=2\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        assert len(lines) == 2
        same_pairwise_elements(
            lines[1].strip().split("\t"),
            "chr21	18842204	chr21	18860085	2	0	++	0,0	18.002	FAIL".split("\t"),
        )


def test_consistent(tmp_path):
    # Run without candidates
    config_file = io.StringIO(f"bam_file={BAM}\noutdir={tmp_path / 'novel'}\nthreads=2\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "novel" / "NAIBR_SVs.bedpe"

    # Generate candidates.bedpe base on results
    new_candidates = tmp_path / "candidates.bedpe"
    with open(new_candidates, "w") as outfile, open(output) as infile:
        line = infile.readlines()[1]
        chr1, pos1, chr2, pos2, _, _, orient, *_ = line.split("\t")
        print(chr1, pos1, pos1, chr2, pos2, pos2, orient, sep="\t", file=outfile)

    # Run with generated candidates
    rerun_config_file = io.StringIO(
        f"bam_file={BAM}\ncandidates={new_candidates}\noutdir={tmp_path / 'rerun'}\nthreads=2\n"
    )
    rerun_configs = Configs.from_file(rerun_config_file)

    exitcode = run(rerun_configs)
    assert exitcode == 0
    output_rerun = tmp_path / "rerun" / "NAIBR_SVs.bedpe"
    with open(output_rerun) as f:
        new_line = f.readlines()[1]
        same_pairwise_elements(line, new_line)
