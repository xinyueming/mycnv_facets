"""Tests for normalizer.py: chromosome statistics logic."""

from __future__ import annotations

import tempfile
from pathlib import Path

from mycnv_facets.normalizer import generate_normal_chrom


class TestGenerateNormalChrom:
    """Tests for normal chromosome identification."""

    def _make_tsv(self, rows: str) -> str:
        """Create a temp TSV file with given data rows."""
        header = "chrom\tstart\tend\talt\tsvtype\tlength\tcopy number\tpanel genes\tt_cn\tl_cn\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write(header + rows)
            f.flush()
            return f.name

    def test_single_chrom_cn2(self):
        """Chromosome appearing once with copy number 2 is normal."""
        tsv = self._make_tsv("chr1\t100\t200\tDUP\tDUP\t100\t2\t.\t2\t1\n")
        out = tsv.replace(".tsv", ".txt")
        result = generate_normal_chrom(tsv, out)

        assert result == "chr1"
        with open(out) as f:
            assert f.read().strip() == "chr1"

    def test_single_chrom_cn3(self):
        """Chromosome with copy number != 2 is not normal."""
        tsv = self._make_tsv("chr1\t100\t200\tDUP\tDUP\t100\t3\t.\t3\t1\n")
        out = tsv.replace(".tsv", ".txt")
        result = generate_normal_chrom(tsv, out)

        assert result == ""

    def test_chrom_appears_twice(self):
        """Chromosome appearing more than once is excluded."""
        tsv = self._make_tsv(
            "chr1\t100\t200\tDUP\tDUP\t100\t2\t.\t2\t1\n"
            "chr1\t300\t400\tDEL\tDEL\t100\t2\t.\t2\t1\n"
        )
        out = tsv.replace(".tsv", ".txt")
        result = generate_normal_chrom(tsv, out)

        assert result == ""

    def test_multiple_normal_chroms(self):
        """Multiple normal chromosomes should be comma-separated."""
        tsv = self._make_tsv(
            "chr1\t100\t200\tDUP\tDUP\t100\t2\t.\t2\t1\n"
            "chr2\t100\t200\tDUP\tDUP\t100\t3\t.\t3\t1\n"
            "chr2\t300\t400\tDEL\tDEL\t100\t3\t.\t3\t0\n"
            "chr3\t100\t200\tDUP\tDUP\t100\t2\t.\t2\t1\n"
        )
        out = tsv.replace(".tsv", ".txt")
        result = generate_normal_chrom(tsv, out)

        # chr1: once, cn=2 → included
        # chr2: twice → excluded
        # chr3: once, cn=2 → included
        assert result == "chr1,chr3"

    def test_missing_tcn(self):
        """Missing copy number ('.') should not be treated as normal."""
        tsv = self._make_tsv("chr1\t100\t200\tDUP\tDUP\t100\t.\t.\t.\t.\n")
        out = tsv.replace(".tsv", ".txt")
        result = generate_normal_chrom(tsv, out)

        assert result == ""

    def test_all_abnormal(self):
        """All chromosomes abnormal should return empty."""
        tsv = self._make_tsv(
            "chr1\t100\t200\tDUP\tDUP\t100\t4\tGENE_A\t4\t1\n"
            "chr2\t100\t200\tcn-LOH\tLOH\t100\t2\t.\t2\t0\n"
            "chr2\t300\t400\tDUP\tDUP\t100\t3\t.\t3\t1\n"
        )
        out = tsv.replace(".tsv", ".txt")
        result = generate_normal_chrom(tsv, out)

        # chr1: once but cn=4 → excluded
        # chr2: appears twice → excluded
        assert result == ""

    def test_output_file_created(self):
        """Output file should always be created."""
        tsv = self._make_tsv("chr1\t100\t200\tDUP\tDUP\t100\t3\t.\t3\t1\n")
        out = tsv.replace(".tsv", ".txt")
        generate_normal_chrom(tsv, out)

        assert Path(out).exists()
