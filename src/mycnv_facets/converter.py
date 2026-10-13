"""Step 3: Convert annotated VCF to a tabular TSV format."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path


def _parse_info_field(info_str: str) -> dict[str, str]:
    """Parse a VCF INFO string into a dict."""
    result: dict[str, str] = {}
    if not info_str or info_str == ".":
        return result
    for token in info_str.split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value
        else:
            result[token] = "1"
    return result


def convert_vcf_to_table(
    vcf_path: str,
    out_tsv_path: str,
) -> str:
    """Convert an annotated VCF to a TSV table.

    Columns: chrom, start, end, alt, svtype, length, copy number, panel genes,
    t_cn, l_cn

    The ``alt`` column uses cn-LOH logic:
    - If TCN_EM == 2 and LCN_EM is not None and TCN_EM - LCN_EM != LCN_EM,
      then alt = cn-LOH
    - Otherwise alt = SVTYPE

    Parameters
    ----------
    vcf_path : str
        Input annotated VCF (plain or gzipped).
    out_tsv_path : str
        Output TSV path.

    Returns
    -------
    str
        Path to the output TSV.
    """
    is_gz_in = str(vcf_path).endswith(".gz")
    is_gz_out = str(out_tsv_path).endswith(".gz")

    if is_gz_in:
        fh_in = gzip.open(vcf_path, "rt")
    else:
        fh_in = open(vcf_path, "r")

    try:
        if is_gz_out:
            tmp_path = out_tsv_path.replace(".gz", "")
            fh_out = open(tmp_path, "w")
        else:
            tmp_path = None
            fh_out = open(out_tsv_path, "w")

        try:
            header = [
                "chrom",
                "start",
                "end",
                "alt",
                "svtype",
                "length",
                "copy number",
                "panel genes",
                "t_cn",
                "l_cn",
            ]
            fh_out.write("\t".join(header) + "\n")

            for line in fh_in:
                if line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                info_str = parts[7]
                info = _parse_info_field(info_str)

                end = int(info.get("END", pos))
                svtype = info.get("SVTYPE", ".")
                svlen_raw = info.get("SVLEN", ".")
                svlen = abs(int(svlen_raw)) if svlen_raw != "." else "."
                tcn_em = info.get("TCN_EM", ".")
                lcn_em = info.get("LCN_EM", ".")

                # Compute alt column
                alt = svtype
                if tcn_em != "." and lcn_em != ".":
                    tcn_val = int(tcn_em)
                    lcn_val = int(lcn_em)
                    if tcn_val == 2 and tcn_val - lcn_val != lcn_val:
                        alt = "cn-LOH"

                copy_number = tcn_em if tcn_em != "." else "."
                panel_genes = info.get("genes_in_panel", ".")

                row = [
                    chrom,
                    str(pos),
                    str(end),
                    alt,
                    svtype,
                    str(svlen),
                    copy_number,
                    panel_genes,
                    tcn_em,
                    lcn_em,
                ]
                fh_out.write("\t".join(row) + "\n")
        finally:
            fh_out.close()

        if tmp_path:
            with open(tmp_path, "rb") as f_in, gzip.open(
                out_tsv_path, "wb"
            ) as f_out:
                shutil.copyfileobj(f_in, f_out)
            Path(tmp_path).unlink()

    finally:
        fh_in.close()

    return out_tsv_path
