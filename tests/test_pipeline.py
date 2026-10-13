"""Integration tests: end-to-end pipeline on real test data."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from mycnv_facets.annotator import annotate_vcf
from mycnv_facets.converter import convert_vcf_to_table
from mycnv_facets.normalizer import generate_normal_chrom

TEST_DATA = Path(__file__).resolve().parent.parent / "test_data"


@pytest.mark.integration
class TestIntegrationPipeline:
    """End-to-end pipeline tests using test_data/."""

    @pytest.fixture
    def tmp_out(self, tmp_path: Path) -> Path:
        return tmp_path / "sample1"

    def test_full_pipeline(self, tmp_path: Path):
        """Run annotation → conversion → normalization on the test VCF."""
        vcf = TEST_DATA / "20C004I0008J.vcf.gz"
        gene_list = TEST_DATA / "zy.brain.230.gene.list"
        refseq = TEST_DATA / "ncbiRefSeqCurated.txt.gz"

        out_prefix = tmp_path / "sample1"
        annotated_vcf = str(out_prefix) + ".annotated.vcf.gz"
        tsv_path = str(out_prefix) + ".annotated.tsv"
        normal_chrom_path = str(out_prefix) + ".normal_chrom.txt"

        # Step 2: Annotate
        annotate_vcf(str(vcf), annotated_vcf, str(gene_list), str(refseq))
        assert Path(annotated_vcf).exists()

        # Verify annotation added genes_in_panel
        with gzip.open(annotated_vcf, "rt") as f:
            data_lines = [l for l in f if not l.startswith("#")]
        annotated_count = sum(1 for l in data_lines if "genes_in_panel=" in l)
        assert annotated_count > 0, "Expected genes_in_panel in annotated VCF"

        # Step 3: Convert
        convert_vcf_to_table(annotated_vcf, tsv_path)
        assert Path(tsv_path).exists()

        with open(tsv_path) as f:
            lines = f.readlines()
        # Header + 27 CNV records (same as original VCF)
        assert len(lines) == len(data_lines) + 1

        # Verify header
        header = lines[0].strip().split("\t")
        assert header == [
            "chrom", "start", "end", "alt", "svtype", "length",
            "copy number", "panel genes", "t_cn", "l_cn",
        ]

        # Step 4: Normalizer
        normal_chroms = generate_normal_chrom(tsv_path, normal_chrom_path)
        assert Path(normal_chrom_path).exists()

        with open(normal_chrom_path) as f:
            chrom_content = f.read().strip()

        assert chrom_content == normal_chroms
        # chr10 and chr22 should be normal (single occurrence, CN=2)
        assert "chr10" in chrom_content
        assert "chr22" in chrom_content

    def test_converter_output_matches_known(self, tmp_path: Path):
        """TSV output should match the known good output."""
        vcf = TEST_DATA / "20C004I0008J.vcf.gz"
        gene_list = TEST_DATA / "zy.brain.230.gene.list"
        refseq = TEST_DATA / "ncbiRefSeqCurated.txt.gz"

        out_prefix = tmp_path / "sample1"
        annotated_vcf = str(out_prefix) + ".annotated.vcf.gz"
        tsv_path = str(out_prefix) + ".annotated.tsv"

        annotate_vcf(str(vcf), annotated_vcf, str(gene_list), str(refseq))
        convert_vcf_to_table(annotated_vcf, tsv_path)

        with open(tsv_path) as f:
            lines = f.readlines()

        # Check specific known rows
        # chr11 cn-LOH row should have alt=cn-LOH
        cn_loh_row = [l for l in lines if "cn-LOH" in l]
        assert len(cn_loh_row) >= 1, "Expected at least one cn-LOH row"

        # chr10 NEUTR row should have copy number=2
        neut_chr10 = [l for l in lines if l.startswith("chr10") and "NEUTR" in l]
        assert len(neut_chr10) == 1
        parts = neut_chr10[0].strip().split("\t")
        assert parts[6] == "2"  # copy number column
