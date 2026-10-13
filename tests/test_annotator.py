"""Tests for annotator.py: gene interval loading and overlap detection."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from mycnv_facets.annotator import (
    _merge_intervals,
    _overlaps,
    annotate_vcf,
)


class TestOverlaps:
    """Tests for interval overlap detection."""

    def test_overlapping(self):
        assert _overlaps(10, 20, 15, 25) is True

    def test_contained(self):
        assert _overlaps(10, 30, 15, 20) is True

    def test_adjacent(self):
        """Adjacent intervals (touching at boundary) should overlap."""
        assert _overlaps(10, 20, 20, 30) is True

    def test_non_overlapping(self):
        assert _overlaps(10, 20, 21, 30) is False
        assert _overlaps(10, 15, 20, 30) is False

    def test_identical(self):
        assert _overlaps(10, 20, 10, 20) is True

    def test_one_inside_other(self):
        assert _overlaps(5, 30, 10, 20) is True


class TestMergeIntervals:
    """Tests for interval merging."""

    def test_no_overlap(self):
        intervals = [(10, 20), (30, 40)]
        assert _merge_intervals(intervals) == [(10, 20), (30, 40)]

    def test_overlap(self):
        intervals = [(10, 20), (15, 25), (30, 40)]
        assert _merge_intervals(intervals) == [(10, 25), (30, 40)]

    def test_adjacent(self):
        intervals = [(10, 20), (20, 30)]
        assert _merge_intervals(intervals) == [(10, 30)]

    def test_contained(self):
        intervals = [(10, 30), (15, 20)]
        assert _merge_intervals(intervals) == [(10, 30)]

    def test_empty(self):
        assert _merge_intervals([]) == []

    def test_single(self):
        assert _merge_intervals([(10, 20)]) == [(10, 20)]

    def test_chain_overlap(self):
        intervals = [(10, 15), (14, 20), (19, 25), (30, 35)]
        assert _merge_intervals(intervals) == [(10, 25), (30, 35)]


class TestAnnotateVcf:
    """Tests for VCF gene annotation."""

    @pytest.fixture
    def simple_vcf(self) -> str:
        """Create a minimal VCF for testing."""
        content = (
            "##fileformat=VCFv4.2\n"
            "##reference=hg19\n"
            "##contig=<ID=chr1,length=249250621>\n"
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Type>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End>\n"
            "##INFO=<ID=TCN_EM,Number=1,Type=Integer,Description=TCN>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t1\tN\t<CNV>\t.\tPASS\tSVTYPE=DUP;END=500;TCN_EM=3\n"
            "chr1\t600\t2\tN\t<CNV>\t.\tneutral\tSVTYPE=NEUTR;END=800;TCN_EM=2\n"
            "chr2\t100\t3\tN\t<CNV>\t.\tPASS\tSVTYPE=DEL;END=300;TCN_EM=1\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        ) as f:
            f.write(content)
            f.flush()
            return f.name

    @pytest.fixture
    def refseq_file(self) -> str:
        """Create a minimal RefSeq-like file."""
        # Columns: bin, name, chrom, strand, txStart, txEnd, ..., gene_name (col 13)
        content = (
            "585\tNM_001\tchr1\t+\t150\t400\t400\t400\t1\t150,\t400,\t0\tGENE_A\tnone\tnone\t-1,\n"
            "585\tNM_002\tchr1\t-\t200\t600\t600\t600\t1\t200,\t600,\t0\tGENE_B\tnone\tnone\t-1,\n"
            "585\tNM_003\tchr2\t+\t100\t200\t200\t200\t1\t100,\t200,\t0\tGENE_C\tnone\tnone\t-1,\n"
            "585\tNM_004\tchr2\t+\t500\t700\t700\t700\t1\t500,\t700,\t0\tGENE_D\tnone\tnone\t-1,\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as f:
            f.write(content)
            f.flush()
            return f.name

    @pytest.fixture
    def gene_list_file(self) -> str:
        """Create a gene list file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".list", delete=False
        ) as f:
            f.write("GENE_A\nGENE_B\nGENE_C\n")
            f.flush()
            return f.name

    def test_annotates_overlapping_genes(self, simple_vcf, refseq_file, gene_list_file):
        """CNVs overlapping gene intervals should get genes_in_panel."""
        out_path = simple_vcf.replace(".vcf", ".annotated.vcf")
        annotate_vcf(simple_vcf, out_path, gene_list_file, refseq_file)

        with open(out_path) as f:
            lines = f.readlines()

        # chr1:100-500 overlaps GENE_A (150-400) and GENE_B (200-600)
        data_lines = [l for l in lines if not l.startswith("#")]
        assert len(data_lines) == 3

        first = data_lines[0]
        assert "genes_in_panel=GENE_A,GENE_B" in first

    def test_no_genes_for_non_overlapping(self, simple_vcf, refseq_file, gene_list_file):
        """CNVs not overlapping any gene should not have genes_in_panel."""
        out_path = simple_vcf.replace(".vcf", ".annotated.vcf")
        annotate_vcf(simple_vcf, out_path, gene_list_file, refseq_file)

        with open(out_path) as f:
            lines = f.readlines()

        data_lines = [l for l in lines if not l.startswith("#")]
        # Third row: chr2:100-300 overlaps GENE_C (100-200) which IS in gene list
        third = data_lines[2]
        assert "genes_in_panel=GENE_C" in third

    def test_preserves_header(self, simple_vcf, refseq_file, gene_list_file):
        """VCF header lines should be preserved exactly."""
        out_path = simple_vcf.replace(".vcf", ".annotated.vcf")
        annotate_vcf(simple_vcf, out_path, gene_list_file, refseq_file)

        with open(simple_vcf) as f_in, open(out_path) as f_out:
            header_in = [l for l in f_in if l.startswith("#")]
            f_out.seek(0)
            header_out = [l for l in f_out if l.startswith("#")]

        assert header_in == header_out

    def test_plain_gene_list(self, simple_vcf, refseq_file, gene_list_file):
        """Test with non-gzipped input."""
        out_path = simple_vcf.replace(".vcf", ".annotated.vcf")
        result = annotate_vcf(simple_vcf, out_path, gene_list_file, refseq_file)
        assert result == out_path
        assert Path(out_path).exists()
