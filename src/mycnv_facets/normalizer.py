"""Step 4: Generate normal_chrom.txt from the annotated TSV."""

from __future__ import annotations

import csv
from pathlib import Path


def generate_normal_chrom(
    tsv_path: str,
    out_path: str,
) -> str:
    """Identify chromosomes with copy number 2 that appear once in the TSV.

    Logic:
    1. Count how many times each chrom appears in the TSV.
    2. For chroms appearing exactly once, check if their ``copy number`` is 2.
    3. Write matching chroms, comma-separated, to the output file.

    Parameters
    ----------
    tsv_path : str
        Path to the annotated TSV (from converter.py).
    out_path : str
        Output file path (normal_chrom.txt).

    Returns
    -------
    str
        The comma-separated normal chromosome string.
    """
    chrom_count: dict[str, int] = {}
    chrom_tcn: dict[str, str] = {}

    opener = open
    is_gz = str(tsv_path).endswith(".gz")
    if is_gz:
        import gzip

        def opener(p: str):
            return gzip.open(p, "rt")

    with opener(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row.get("chrom", "").strip()
            tcn = row.get("copy number", ".").strip()
            if not chrom:
                continue
            chrom_count[chrom] = chrom_count.get(chrom, 0) + 1
            chrom_tcn[chrom] = tcn

    normal_chroms = []
    for chrom, count in chrom_count.items():
        if count == 1 and chrom_tcn.get(chrom) == "2":
            normal_chroms.append(chrom)

    result = ",".join(normal_chroms)
    Path(out_path).write_text(result + "\n" if result else "\n")
    return result
