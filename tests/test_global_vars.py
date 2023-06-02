from io import StringIO
from pathlib import Path

import pytest

from naibr.global_vars import Configs

BAM = Path("tests/data/example_chr21.bam").absolute()


def test_configs_from_file():
    d = 5000
    mapq = 20
    k = 3
    sd_mult = 1
    min_sv = 900
    min_len = 1000
    max_len = 100_000
    min_reads = 0
    min_discs = 4
    s = "# Some comment\n"  # Skipped
    s += f"bam_file={BAM}\n"
    s += f"d={d}\n"
    s += f"k={k}\n"
    s += f"sd_mult={sd_mult}\n"
    s += f"min_sv={min_sv}\n"
    s += f"min_mapq={mapq}\n"
    s += f"min_len={min_len}\n"
    s += f"max_len={max_len}\n"
    s += f"min_reads={min_reads}\n"
    s += f"min_discs={min_discs}\n"
    s += "outdir=results\n"

    with StringIO(s) as f:
        configs = Configs.from_file(f)

    assert BAM.samefile(configs.BAM_FILE)
    assert configs.MAX_LINKED_DIST == d
    assert configs.MIN_MAPQ == mapq
    assert configs.MIN_BC_OVERLAP == k
    assert configs.SD_MULT == sd_mult
    assert configs.MIN_SV == min_sv
    assert configs.MIN_LEN == min_len
    assert configs.MAX_LEN == max_len
    assert configs.MIN_READS == min_reads
    assert configs.MIN_DISCS == min_discs


def test_missing_bam_file_raises_error(tmp_path):
    config_file_without_bam = StringIO("")
    with pytest.raises(SystemExit):
        _ = Configs.from_file(config_file_without_bam)


def test_empty_prefix_is_set_to_bam_file_name(tmp_path):
    config_file_without_prefix = StringIO(f"bam_file={BAM}\nprefix=\n")
    configs = Configs.from_file(config_file_without_prefix)
    assert configs.PREFIX == BAM.stem


def test_set_prefix_in_config(tmp_path):
    prefix = "my_prefix"
    config_file_with_prefix = StringIO(f"bam_file={BAM}\nprefix={prefix}\n")
    configs = Configs.from_file(config_file_with_prefix)
    assert configs.PREFIX == prefix
