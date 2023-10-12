from io import StringIO

import pytest

from naibr.utils import input_candidates, UnionDict


def test_input_candates_bepde_format():
    input_cands = [("chr1", 1000, "chr1", 4000, "+-"), ("chr1", 5000, "chr1", 13000, "+-")]
    # Convert to BEDPE format
    input = "\n".join([f"{c[0]}\t{c[1]}\t{c[1]}\t{c[2]}\t{c[3]}\t{c[3]}\t{c[4]}" for c in input_cands])
    with StringIO(input) as f:
        cands = input_candidates(f)
        assert len(cands) == 2
        cand1, cand2 = cands
        assert cand1 == input_cands[0]
        assert cand2 == input_cands[1]


def test_input_candates_naibr_format():
    input_cands = [("chr1", 1000, "chr1", 4000, "+-"), ("chr1", 5000, "chr1", 13000, "+-")]
    # Convert to NAIBR format with header
    input = (
        "Chr1\tBreak1\tChr2\tBreak2\tSplit molecules\tDiscordant reads\tOrientation\tHaplotype\tScore\t"
        "Pass filter\n"
    )
    input += "\n".join(
        [f"{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\t30\t6\t{c[4]}\t2,2\t345.980\tPASS" for c in input_cands]
    )
    print(input)
    with StringIO(input) as f:
        cands = input_candidates(f)
        assert len(cands) == 2
        cand1, cand2 = cands
        assert cand1 == input_cands[0]
        assert cand2 == input_cands[1]


def test_input_candates_filter_too_short():
    s = "chr1\t1000\t1000\tchr1\t1100\t1100\t--\n"  # too short
    s += "chr1\t2000\t2000\tchr1\t4000\t4000\t+-\n"
    s += "chr1\t10000\t10000\tchr1\t14000\t14000\t-+\n"
    with StringIO(s) as f:
        assert len(input_candidates(f, min_sv=1000)) == 2


def test_union_dict_combine_lists():
    d1 = UnionDict(list)
    d1["a"].extend([1, 2, 3])
    d1["b"].extend([7, 8, 9])

    d2 = UnionDict(list)
    d2["a"].extend([3, 4, 5])
    d2["b"].extend([8, 9, 10])

    d1.combine(d2)
    assert sorted(d1["a"]) == sorted([1, 2, 3, 3, 4, 5])
    assert sorted(d1["b"]) == sorted([7, 8, 8, 9, 9, 10])


def test_union_dict_combine_sets():
    d1 = UnionDict(set)
    d1["a"].union({1, 2, 3})
    d1["b"].union({7, 8, 9})

    d2 = UnionDict(set)
    d2["a"].union({3, 4, 5})
    d2["b"].union({8, 9, 10})

    d1.combine(d2)
    assert d1["a"] - {1, 2, 3, 4, 5} == set()
    assert d1["b"] - {7, 8, 9, 10} == set()


def test_union_dict_combine_different_fails():
    d1 = UnionDict(set)
    d1["a"].union({1, 2, 3})
    d1["b"].union({7, 8, 9})

    d2 = UnionDict(list)
    d2["a"].extend([3, 4, 5])
    d2["b"].extend([8, 9, 10])

    with pytest.raises(ValueError):
        d1.combine(d2)

    with pytest.raises(ValueError):
        d2.combine(d1)
