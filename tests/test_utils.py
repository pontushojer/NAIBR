from io import StringIO

from naibr.utils import input_candidates


def test_input_candates_filter_too_short():
    s = "chr1\t1000\t1000\tchr1\t1100\t1100\t--\n"  # too short
    s += "chr1\t2000\t2000\tchr1\t4000\t4000\t+-\n"
    s += "chr1\t10000\t10000\tchr1\t14000\t14000\t-+\n"
    with StringIO(s) as f:
        assert len(input_candidates(f, min_sv=1000)) == 2
