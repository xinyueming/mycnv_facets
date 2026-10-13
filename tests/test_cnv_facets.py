"""Tests for cnv_facets.py: parameter conversion and command building."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from mycnv_facets.cnv_facets import find_cnv_facets_r, run_cnv_facets


class TestFindCnvFacetsR:
    """Tests for locating the cnv_facets.R script."""

    def test_raises_when_not_found(self):
        with patch("mycnv_facets.cnv_facets.Path.exists", return_value=False), \
             patch("mycnv_facets.cnv_facets.shutil.which", return_value=None), \
             pytest.raises(FileNotFoundError, match="cnv_facets.R not found"):
            find_cnv_facets_r()


class TestRunCnvFacets:
    """Tests for the run_cnv_facets command builder."""

    def test_basic_command(self):
        """Test that basic required arguments are passed correctly."""
        with patch("mycnv_facets.cnv_facets.find_cnv_facets_r", return_value="/fake/cnv_facets.R"), \
             patch("mycnv_facets.cnv_facets.Path.exists", return_value=False), \
             patch("mycnv_facets.cnv_facets.shutil.which", return_value="/usr/bin/Rscript"), \
             patch("subprocess.run") as mock_run:

            mock_run.return_value = MagicMock(returncode=0)
            run_cnv_facets(
                tumour_bam="t.bam",
                normal_bam="n.bam",
                snp_vcf="snps.vcf.gz",
                out="/tmp/sample1",
            )

            cmd = mock_run.call_args[0][0]
            assert "Rscript" in cmd[0] or cmd[0] == "Rscript"
            assert "-o" in cmd
            assert "/tmp/sample1" in cmd
            assert "-t" in cmd
            assert "t.bam" in cmd
            assert "-n" in cmd
            assert "n.bam" in cmd
            assert "-vcf" in cmd
            assert "snps.vcf.gz" in cmd

    def test_default_values(self):
        """Test default argument values."""
        with patch("mycnv_facets.cnv_facets.find_cnv_facets_r", return_value="/fake/cnv_facets.R"), \
             patch("subprocess.run") as mock_run:

            mock_run.return_value = MagicMock(returncode=0)
            run_cnv_facets(
                tumour_bam="t.bam",
                normal_bam="n.bam",
                snp_vcf="snps.vcf.gz",
                out="/tmp/sample1",
            )

            cmd = mock_run.call_args[0][0]
            # Default mapq=5
            assert "-mq" in cmd
            assert "5" in cmd
            # Default baq=10
            assert "-bq" in cmd
            assert "10" in cmd
            # Default depth=(25, 4000)
            assert "25" in cmd
            assert "4000" in cmd
            # Default cval=(25, 150)
            assert "150" in cmd
            # Default nbhd_snp="auto"
            assert "auto" in cmd
            # Default gbuild="hg19"
            assert "hg19" in cmd

    def test_custom_values(self):
        """Test custom argument values are passed through."""
        with patch("mycnv_facets.cnv_facets.find_cnv_facets_r", return_value="/fake/cnv_facets.R"), \
             patch("subprocess.run") as mock_run:

            mock_run.return_value = MagicMock(returncode=0)
            run_cnv_facets(
                tumour_bam="t.bam",
                normal_bam="n.bam",
                snp_vcf="snps.vcf.gz",
                out="/tmp/sample1",
                mapq=25,
                baq=20,
                depth=(100, 2000),
                cval=(50, 200),
                nbhd_snp="300",
                gbuild="hg38",
                targets="panel.bed",
                annotation="genes.bed",
                unmatched=True,
                no_cov_plot=True,
                count_orphans=True,
                nprocs=4,
            )

            cmd = mock_run.call_args[0][0]
            cmd_str = " ".join(cmd)
            assert "25" in cmd_str
            assert "20" in cmd_str
            assert "100" in cmd_str
            assert "2000" in cmd_str
            assert "hg38" in cmd_str
            assert "panel.bed" in cmd_str
            assert "genes.bed" in cmd_str
            assert "-A" in cmd
            assert "-u" in cmd
            assert "-np" in cmd

    def test_returns_vcf_path(self):
        """Test that the function returns the expected VCF path."""
        with patch("mycnv_facets.cnv_facets.find_cnv_facets_r", return_value="/fake/cnv_facets.R"), \
             patch("subprocess.run") as mock_run:

            mock_run.return_value = MagicMock(returncode=0)
            result = run_cnv_facets(
                tumour_bam="t.bam",
                normal_bam="n.bam",
                snp_vcf="snps.vcf.gz",
                out="/tmp/sample1",
            )
            assert result == "/tmp/sample1.vcf.gz"

    def test_pileup_mode(self):
        """Test pileup input instead of BAM files."""
        with patch("mycnv_facets.cnv_facets.find_cnv_facets_r", return_value="/fake/cnv_facets.R"), \
             patch("subprocess.run") as mock_run:

            mock_run.return_value = MagicMock(returncode=0)
            run_cnv_facets(
                pileup="existing.csv.gz",
                out="/tmp/sample1",
            )

            cmd = mock_run.call_args[0][0]
            assert "-p" in cmd
            assert "existing.csv.gz" in cmd
