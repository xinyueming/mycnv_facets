"""Tests for mycnv_facets core modules."""

from __future__ import annotations

import os
import tempfile

from mycnv_facets.annotator import _merge_intervals, _overlaps
from mycnv_facets.converter import _parse_info_field
from mycnv_facets.normalizer import generate_normal_chrom


def test_parse_info_field():
    info = "SVTYPE=DUP;SVLEN=100;END=200;TCN_EM=3;LCN_EM=1"
    result = _parse_info_field(info)
    assert result["SVTYPE"] == "DUP"
    assert result["SVLEN"] == "100"
    assert result["END"] == "200"
    assert result["TCN_EM"] == "3"
    assert result["LCN_EM"] == "1"


def test_parse_info_empty():
    assert _parse_info_field(".") == {}
    assert _parse_info_field("") == {}


def test_overlaps():
    assert _overlaps(10, 20, 15, 25) is True
    assert _overlaps(10, 20, 21, 30) is False
    assert _overlaps(10, 20, 10, 20) is True
    assert _overlaps(10, 20, 5, 15) is True


def test_merge_intervals():
    intervals = [(10, 20), (15, 25), (30, 40)]
    merged = _merge_intervals(intervals)
    assert merged == [(10, 25), (30, 40)]


def test_merge_intervals_adjacent():
    intervals = [(10, 20), (20, 30)]
    merged = _merge_intervals(intervals)
    assert merged == [(10, 30)]


def test_generate_normal_chrom():
    """Test normal_chrom generation with synthetic TSV data."""
    tsv_content = (
        "chrom\tstart\tend\talt\tsvtype\tlength\tcopy number\tpanel genes\tt_cn\tl_cn\n"
        "chr1\t100\t200\tDUP\tDUP\t100\t2\t.\t2\t1\n"
        "chr2\t100\t200\tDUP\tDUP\t100\t3\t.\t3\t1\n"
        "chr2\t300\t400\tDEL\tDEL\t100\t3\t.\t3\t1\n"
        "chr3\t100\t200\tDUP\tDUP\t100\t3\t.\t3\t1\n"
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv_content)
        f.flush()
        tsv_path = f.name

    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            out_path = f.name

        result = generate_normal_chrom(tsv_path, out_path)
        with open(out_path) as f:
            content = f.read().strip()

        # chr1 appears once with copy number 2 -> included
        # chr2 appears twice -> excluded
        # chr3 appears once with copy number 3 -> excluded
        assert result == "chr1"
        assert content == "chr1"
        os.unlink(out_path)
    finally:
        os.unlink(tsv_path)
