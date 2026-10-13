"""Step 2: Annotate VCF CNV intervals with gene panel membership."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or adjacent intervals."""
    if not intervals:
        return []
    sorted_iv = sorted(intervals)
    merged = [list(sorted_iv[0])]
    for start, end in sorted_iv[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def _overlaps(s1: int, e1: int, s2: int, e2: int) -> bool:
    """Check if two intervals overlap."""
    return max(s1, s2) <= min(e1, e2)


def annotate_vcf(
    vcf_path: str,
    out_vcf_path: str,
    gene_list_path: str,
    refseq_path: str,
) -> str:
    """Annotate a VCF with ``genes_in_panel`` INFO field.

    Reads the input VCF (plain or gzipped), looks up each CNV interval against
    the gene panel, and writes an annotated VCF.

    Parameters
    ----------
    vcf_path : str
        Input VCF (can be .vcf or .vcf.gz).
    out_vcf_path : str
        Output annotated VCF.
    gene_list_path : str
        One gene symbol per line.
    refseq_path : str
        ncbiRefSeqCurated formatted file.

    Returns
    -------
    str
        Path to the output VCF.
    """
    # Load gene list
    with open(gene_list_path) as fh:
        gene_set = {line.strip() for line in fh if line.strip()}

    # Build per-gene interval index: gene -> [(chrom, start, end), ...]
    gene_intervals: dict[str, list[tuple[str, int, int]]] = {}

    opener = gzip.open if str(refseq_path).endswith(".gz") else open
    with opener(refseq_path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13:
                continue
            chrom = parts[2]
            tx_start = int(parts[4])
            tx_end = int(parts[5])
            gene_name = parts[12]
            if gene_name not in gene_set:
                continue
            if gene_name not in gene_intervals:
                gene_intervals[gene_name] = []
            gene_intervals[gene_name].append((chrom, tx_start, tx_end))

    # Merge per-gene intervals per chrom for faster lookup
    # Then build a flat chrom -> [(start, end, gene_name), ...] index
    chrom_index: dict[str, list[tuple[int, int, str]]] = {}
    for gene, intervals in gene_intervals.items():
        # Group by chrom
        by_chrom: dict[str, list[tuple[int, int]]] = {}
        for chrom, s, e in intervals:
            by_chrom.setdefault(chrom, []).append((s, e))
        for chrom, ivs in by_chrom.items():
            merged = _merge_intervals(ivs)
            chrom_index.setdefault(chrom, []).extend(
                [(s, e, gene) for s, e in merged]
            )

    # Open input
    is_gz_in = str(vcf_path).endswith(".gz")
    is_gz_out = str(out_vcf_path).endswith(".gz")

    if is_gz_in:
        fh_in = gzip.open(vcf_path, "rt")
    else:
        fh_in = open(vcf_path, "r")

    try:
        if is_gz_out:
            # Write to temp plain file, then gzip
            tmp_path = out_vcf_path.replace(".gz", "")
            fh_out = open(tmp_path, "w")
        else:
            tmp_path = None
            fh_out = open(out_vcf_path, "w")

        try:
            for line in fh_in:
                if line.startswith("#"):
                    fh_out.write(line)
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    fh_out.write(line)
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                info = parts[7]

                # Parse END from INFO
                end = pos
                for token in info.split(";"):
                    if token.startswith("END="):
                        end = int(token.split("=", 1)[1])
                        break

                # Find overlapping genes
                genes = []
                if chrom in chrom_index:
                    for gs, ge, gname in chrom_index[chrom]:
                        if _overlaps(pos, end, gs, ge):
                            genes.append(gname)

                # Append genes_in_panel to INFO
                if genes:
                    # Escape commas and special chars in gene names for VCF
                    gene_str = ",".join(genes)
                    if info == "." or info == "":
                        new_info = f"genes_in_panel={gene_str}"
                    else:
                        new_info = f"{info};genes_in_panel={gene_str}"
                    parts[7] = new_info
                    fh_out.write("\t".join(parts) + "\n")
                else:
                    fh_out.write(line)
        finally:
            fh_out.close()

        if tmp_path:
            with open(tmp_path, "rb") as f_in, gzip.open(
                out_vcf_path, "wb"
            ) as f_out:
                shutil.copyfileobj(f_in, f_out)
            Path(tmp_path).unlink()

    finally:
        fh_in.close()

    return out_vcf_path
