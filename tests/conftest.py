"""Shared fixtures and test data paths."""

from __future__ import annotations

from pathlib import Path

import pytest

TEST_DATA = Path(__file__).resolve().parent.parent / "test_data"


@pytest.fixture
def vcf_path() -> Path:
    return TEST_DATA / "20C004I0008J.vcf.gz"


@pytest.fixture
def gene_list_path() -> Path:
    return TEST_DATA / "zy.brain.230.gene.list"


@pytest.fixture
def refseq_path() -> Path:
    return TEST_DATA / "ncbiRefSeqCurated.txt.gz"
