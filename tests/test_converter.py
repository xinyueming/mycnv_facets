"""Tests for converter.py: VCF parsing and table generation."""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

import pytest

from mycnv_facets.converter import _parse_info_field, convert_vcf_to_table


class TestParseInfoField:
    """Tests for VCF INFO field parsing."""

    def test_single_key_value(self):
        result = _parse_info_field("SVTYPE=DUP")
        assert result == {"SVTYPE": "DUP"}

    def test_multiple_key_values(self):
        result = _parse_info_field("SVTYPE=DUP;END=500;TCN_EM=3;LCN_EM=1")
        assert result["SVTYPE"] == "DUP"
        assert result["END"] == "500"
        assert result["TCN_EM"] == "3"
        assert result["LCN_EM"] == "1"

    def test_empty(self):
        assert _parse_info_field(".") == {}
        assert _parse_info_field("") == {}

    def test_flag_without_value(self):
        result = _parse_info_field("SOME_FLAG")
        assert result["SOME_FLAG"] == "1"

    def test_mixed(self):
        result = _parse_info_field("SVTYPE=DEL;SOME_FLAG;END=100")
        assert result["SVTYPE"] == "DEL"
        assert result["SOME_FLAG"] == "1"
        assert result["END"] == "100"


class TestConvertVcfToTable:
    """Tests for VCF to TSV conversion."""

    @pytest.fixture
    def sample_vcf(self) -> str:
        """Create a minimal annotated VCF for testing."""
        content = (
            "##fileformat=VCFv4.2\n"
            "##reference=hg19\n"
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Type>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End>\n"
            "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=Length>\n"
            "##INFO=<ID=TCN_EM,Number=1,Type=Integer,Description=TCN>\n"
            "##INFO=<ID=LCN_EM,Number=1,Type=Integer,Description=LCN>\n"
            "##INFO=<ID=genes_in_panel,Number=.,Type=String,Description=Genes>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t1\tN\t<CNV>\t.\tPASS\tSVTYPE=DUP;SVLEN=400;END=500;TCN_EM=4;LCN_EM=1;genes_in_panel=GENE_A,GENE_B\n"
            "chr2\t200\t2\tN\t<CNV>\t.\tneutral\tSVTYPE=NEUTR;SVLEN=100;END=300;TCN_EM=2;LCN_EM=1\n"
            "chr11\t300\t3\tN\t<CNV>\t.\tPASS\tSVTYPE=LOH;SVLEN=100;END=400;TCN_EM=2;LCN_EM=0;genes_in_panel=HRAS\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        ) as f:
            f.write(content)
            f.flush()
            return f.name

    def test_column_headers(self, sample_vcf):
        """TSV output should have expected column headers."""
        out_path = sample_vcf.replace(".vcf", ".tsv")
        convert_vcf_to_table(sample_vcf, out_path)

        with open(out_path) as f:
            header = f.readline().strip().split("\t")

        expected = [
            "chrom", "start", "end", "alt", "svtype", "length",
            "copy number", "panel genes", "t_cn", "l_cn",
        ]
        assert header == expected

    def test_row_count(self, sample_vcf):
        """TSV should have one row per VCF data record."""
        out_path = sample_vcf.replace(".vcf", ".tsv")
        convert_vcf_to_table(sample_vcf, out_path)

        with open(out_path) as f:
            lines = [l for l in f if l.strip()]

        assert len(lines) == 4  # 1 header + 3 data rows

    def test_field_mapping(self, sample_vcf):
        """VCF fields should be mapped to correct TSV columns."""
        out_path = sample_vcf.replace(".vcf", ".tsv")
        convert_vcf_to_table(sample_vcf, out_path)

        with open(out_path) as f:
            lines = f.readlines()

        # First data row: chr1, DUP, copy number=4
        row1 = lines[1].strip().split("\t")
        assert row1[0] == "chr1"    # chrom
        assert row1[1] == "100"     # start
        assert row1[2] == "500"     # end
        assert row1[3] == "DUP"     # alt
        assert row1[4] == "DUP"     # svtype
        assert row1[5] == "400"     # length
        assert row1[6] == "4"       # copy number
        assert row1[7] == "GENE_A,GENE_B"  # panel genes
        assert row1[8] == "4"       # t_cn
        assert row1[9] == "1"       # l_cn

    def test_cn_loh_logic(self, sample_vcf):
        """TCN_EM=2, LCN_EM=0 should produce alt=cn-LOH."""
        out_path = sample_vcf.replace(".vcf", ".tsv")
        convert_vcf_to_table(sample_vcf, out_path)

        with open(out_path) as f:
            lines = f.readlines()

        # Third data row: chr11, LOH, TCN_EM=2, LCN_EM=0
        # TCN_EM - LCN_EM = 2 - 0 = 2 != LCN_EM (0) => cn-LOH
        row3 = lines[3].strip().split("\t")
        assert row3[3] == "cn-LOH"   # alt should be cn-LOH
        assert row3[4] == "LOH"      # svtype remains LOH

    def test_neutr_not_cn_loh(self, sample_vcf):
        """TCN_EM=2, LCN_EM=1 should NOT produce cn-LOH."""
        out_path = sample_vcf.replace(".vcf", ".tsv")
        convert_vcf_to_table(sample_vcf, out_path)

        with open(out_path) as f:
            lines = f.readlines()

        # Second data row: chr2, NEUTR, TCN_EM=2, LCN_EM=1
        # TCN_EM - LCN_EM = 2 - 1 = 1 == LCN_EM => not cn-LOH
        row2 = lines[2].strip().split("\t")
        assert row2[3] == "NEUTR"   # alt remains NEUTR

    def test_tcn_em_less_than_2_alt_is_del(self):
        """TCN_EM < 2 should produce alt=DEL regardless of SVTYPE."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t1\tN\t<CNV>\t.\tPASS\tSVTYPE=DUP;SVLEN=100;END=200;TCN_EM=1;LCN_EM=0\n"
            "chr2\t200\t2\tN\t<CNV>\t.\tPASS\tSVTYPE=HEMIZYG;SVLEN=50;END=250;TCN_EM=0;LCN_EM=0\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(content)
            f.flush()
            vcf_path = f.name

        out_path = vcf_path.replace(".vcf", ".tsv")
        convert_vcf_to_table(vcf_path, out_path)

        with open(out_path) as f:
            lines = f.readlines()

        # TCN_EM=1 with SVTYPE=DUP → alt=DEL
        row1 = lines[1].strip().split("\t")
        assert row1[3] == "DEL"
        assert row1[4] == "DUP"   # svtype column still shows original

        # TCN_EM=0 with SVTYPE=HEMIZYG → alt=DEL
        row2 = lines[2].strip().split("\t")
        assert row2[3] == "DEL"
        assert row2[4] == "HEMIZYG"

    def test_missing_panel_genes(self, sample_vcf):
        """Records without genes_in_panel should show '.'."""
        out_path = sample_vcf.replace(".vcf", ".tsv")
        convert_vcf_to_table(sample_vcf, out_path)

        with open(out_path) as f:
            lines = f.readlines()

        # Second row has no genes_in_panel
        row2 = lines[2].strip().split("\t")
        assert row2[7] == "."

    def test_gzipped_output(self, sample_vcf):
        """Output should be gzipped when path ends with .gz."""
        out_path = sample_vcf.replace(".vcf", ".tsv.gz")
        convert_vcf_to_table(sample_vcf, out_path)

        assert Path(out_path).exists()
        with gzip.open(out_path, "rt") as f:
            header = f.readline().strip().split("\t")
        expected = ["chrom", "start", "end", "alt", "svtype", "length",
                    "copy number", "panel genes", "t_cn", "l_cn"]
        assert header == expected
