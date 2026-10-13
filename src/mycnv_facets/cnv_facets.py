"""Step 1: Run cnv_facets.R to detect CNVs."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def find_cnv_facets_r() -> str:
    """Locate the cnv_facets.R script.

    Checks in order:
    1. The bin/ directory next to this package
    2. On PATH via ``which cnv_facets.R``
    """
    # Try bin/ relative to this package
    here = Path(__file__).resolve().parent
    candidate = here.parent.parent / "bin" / "cnv_facets.R"
    if candidate.exists():
        return str(candidate)

    # Try on PATH
    path = shutil.which("cnv_facets.R")
    if path:
        return path

    raise FileNotFoundError(
        "cnv_facets.R not found. Install it or place bin/cnv_facets.R next to this package."
    )


def run_cnv_facets(
    *,
    tumour_bam: str | None = None,
    normal_bam: str | None = None,
    snp_vcf: str | None = None,
    out: str,
    pileup: str | None = None,
    mapq: int = 5,
    baq: int = 10,
    count_orphans: bool = False,
    nprocs: int = 1,
    depth: tuple[int, int] = (25, 4000),
    targets: str | None = None,
    cval: tuple[float, float] = (25, 150),
    nbhd_snp: str = "auto",
    annotation: str | None = None,
    gbuild: str = "hg19",
    unmatched: bool = False,
    no_cov_plot: bool = False,
    rnd_seed: str | None = None,
) -> str:
    """Run cnv_facets.R and return the output VCF path.

    Parameters match the R script's arguments.  The returned path is
    ``<out>.vcf.gz``.
    """
    rscript = find_cnv_facets_r()

    cmd = ["Rscript", rscript, "-o", out]

    if tumour_bam:
        cmd += ["-t", tumour_bam]
    if normal_bam:
        cmd += ["-n", normal_bam]
    if snp_vcf:
        cmd += ["-vcf", snp_vcf]
    if pileup:
        cmd += ["-p", pileup]

    cmd += ["-mq", str(mapq)]
    cmd += ["-bq", str(baq)]
    if count_orphans:
        cmd.append("-A")
    cmd += ["-N", str(nprocs)]
    cmd += ["-d", str(depth[0]), str(depth[1])]

    if targets:
        cmd += ["-T", targets]
    cmd += ["-cv", str(cval[0]), str(cval[1])]
    cmd += ["-snp", nbhd_snp]

    if annotation:
        cmd += ["-a", annotation]
    cmd += ["-g", gbuild]
    if unmatched:
        cmd.append("-u")
    if no_cov_plot:
        cmd.append("-np")
    if rnd_seed:
        cmd += ["-s", rnd_seed]

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"cnv_facets.R failed with exit code {result.returncode}:\n{result.stderr}"
        )
    return f"{out}.vcf.gz"
