import io
from pathlib import Path

import pytest

from naibr.__main__ import main, run
from naibr.global_vars import Configs
from naibr.utils import NovelAdjacency

BAM_INV = Path("tests/data/HCC1954T_10xGenomics_chr21_INV.bam").absolute()
CANDIDATES_INV = Path("tests/data/HCC1954T_10xGenomics_chr21_INV.candidates.bedpe").absolute()
CANDIDATES_REV_INV = Path("tests/data/HCC1954T_10xGenomics_chr21_INV.candidates_reversed.bedpe").absolute()
BAM_DEL = Path("tests/data/NA24385_10xGenomics_chr1_DEL.bam").absolute()
BAM_DUP = Path("tests/data/HCC1954T_10xGenomics_chr3_DUP.bam").absolute()
BAM_TRA = Path("tests/data/HCC1954T_10xGenomics_chr_12_20_TRA.bam").absolute()
BAM_SMALL = Path("tests/data/HCC1954T_10xGenomics_chr21_INV.one_barcode.bam").absolute()


def same_pairwise_elements(list1, list2):
    assert len(list1) == len(list2)
    for l1, l2 in zip(list1, list2):
        assert l1 == l2


def iter_novel_adjacencies(bedpe_reader):
    for entry in bedpe_reader:
        # Skip header
        if entry.startswith("Chr1"):
            continue

        yield NovelAdjacency.from_naibr(entry)


def test_help_return_zero():
    for arg in ["-h", "--help", "help"]:
        with pytest.raises(SystemExit):
            errorcode = main([arg])
            assert errorcode == 0


def test_version_return_zero():
    for arg in ["version", "-v", "--version"]:
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


def test_run_from_file(tmp_path):
    configfile = tmp_path / "naibr.config"
    outdir = tmp_path / "results"
    with open(configfile, "w") as f:
        print(f"bam_file={BAM_INV}\ncandidates={CANDIDATES_INV}\noutdir={outdir}\n", file=f)

    errorcode = main([configfile])
    assert errorcode == 0
    assert outdir.exists() and outdir.is_dir()
    assert Path(outdir / "NAIBR_SVs.bedpe").exists()
    assert Path(outdir / "NAIBR_SVs.reformat.bedpe").exists()
    assert Path(outdir / "NAIBR_SVs.vcf").exists()
    assert Path(outdir / "NAIBR_SVs.log").exists()


@pytest.mark.parametrize("threads", [1, 2], ids=["singlethread", "multithread"])
def test_candidates(tmp_path, threads):
    config_file = io.StringIO(
        f"bam_file={BAM_INV}\ncandidates={CANDIDATES_INV}\noutdir={tmp_path}\nthreads={threads}\nDEBUG=True\n"
    )
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as bedpe_reader:
        nas = list(iter_novel_adjacencies(bedpe_reader))
        assert len(nas) == 1
        na = nas[0]
        assert na.pass_threshold
        assert na.orient == "--"
        assert na.chrm1 == na.chrm2 == "chr21"
        assert na.break1 == 18759465
        assert na.break2 == 18906782


def test_blacklist(tmp_path):
    # Create a blacklist region that overlaps with the candidate region
    bed = tmp_path / "blacklist.bed"
    with open(bed, "w") as f:
        print("chr21\t18000000\t19000000", file=f)

    config_file = io.StringIO(f"bam_file={BAM_INV}\nblacklist={bed}\noutdir={tmp_path}\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as f:
        lines = f.readlines()

        # Only header line should be present
        assert len(lines) == 1


def test_candidates_reversed(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM_INV}\ncandidates={CANDIDATES_REV_INV}\noutdir={tmp_path}\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as bedpe_reader:
        nas = list(iter_novel_adjacencies(bedpe_reader))
        assert len(nas) == 1
        na = nas[0]
        assert na.pass_threshold
        assert na.orient == "--"
        assert na.chrm1 == na.chrm2 == "chr21"
        assert na.break1 == 18759465
        assert na.break2 == 18906782


@pytest.mark.parametrize("threads", [1, 2], ids=["singlethread", "multithread"])
def test_call_inversion(tmp_path, threads):
    config_file = io.StringIO(f"bam_file={BAM_INV}\noutdir={tmp_path}\nthreads={threads}\nDEBUG=True\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as bedpe_reader:
        nas = list(iter_novel_adjacencies(bedpe_reader))
        assert len(nas) > 0
        nas_pass = [na for na in nas if na.pass_threshold]
        assert len(nas_pass) == 1
        assert nas_pass[0].orient == "--"
        assert nas_pass[0].chrm1 == nas_pass[0].chrm2 == "chr21"
        assert abs(nas_pass[0].break1 - 18759465) < 500
        assert abs(nas_pass[0].break2 - 18906782) < 500


def test_consistent(tmp_path):
    # Run without candidates
    config_file = io.StringIO(f"bam_file={BAM_INV}\noutdir={tmp_path / 'novel'}\nthreads=2\n")
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
        f"bam_file={BAM_INV}\ncandidates={new_candidates}\noutdir={tmp_path / 'rerun'}\nthreads=2\n"
    )
    rerun_configs = Configs.from_file(rerun_config_file)

    exitcode = run(rerun_configs)
    assert exitcode == 0
    output_rerun = tmp_path / "rerun" / "NAIBR_SVs.bedpe"
    with open(output_rerun) as f:
        new_line = f.readlines()[1]
        same_pairwise_elements(line, new_line)


def test_call_deletion(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM_DEL}\noutdir={tmp_path}\nDEBUG=True\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as bedpe_reader:
        nas = list(iter_novel_adjacencies(bedpe_reader))
        assert len(nas) > 0
        nas_pass = [na for na in nas if na.pass_threshold]
        assert len(nas_pass) == 1
        assert nas_pass[0].orient == "+-"
        assert nas_pass[0].chrm1 == nas_pass[0].chrm2 == "chr1"
        assert abs(nas_pass[0].break1 - 115686862) < 500
        assert abs(nas_pass[0].break2 - 115690219) < 500


def test_call_duplication(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM_DUP}\noutdir={tmp_path}\nDEBUG=True\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as bedpe_reader:
        nas = list(iter_novel_adjacencies(bedpe_reader))
        assert len(nas) > 0
        nas_pass = [na for na in nas if na.pass_threshold]
        assert len(nas_pass) == 1
        assert nas_pass[0].orient == "-+"
        assert nas_pass[0].chrm1 == nas_pass[0].chrm2 == "chr3"
        assert abs(nas_pass[0].break1 - 101569210) < 500
        assert abs(nas_pass[0].break2 - 101663759) < 500


def test_call_interchromosomal(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM_TRA}\noutdir={tmp_path}\nDEBUG=True\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as bedpe_reader:
        nas = list(iter_novel_adjacencies(bedpe_reader))
        assert len(nas) > 0
        nas_pass = [na for na in nas if na.pass_threshold]
        assert len(nas_pass) == 1
        assert nas_pass[0].orient == "--"
        assert nas_pass[0].chrm1 == "chr12"
        assert nas_pass[0].chrm2 == "chr20"
        assert abs(nas_pass[0].break1 - 40393752) < 500
        assert abs(nas_pass[0].break2 - 53974890) < 500


def test_too_few_barcodes(tmp_path):
    config_file = io.StringIO(f"bam_file={BAM_SMALL}\noutdir={tmp_path}\nDEBUG=True\n")
    configs = Configs.from_file(config_file)

    exitcode = run(configs)
    assert exitcode == 0
    output = tmp_path / "NAIBR_SVs.bedpe"
    with open(output) as bedpe_reader:
        nas = list(iter_novel_adjacencies(bedpe_reader))
        assert len(nas) == 0
