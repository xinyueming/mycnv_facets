"""Command-line interface for mycnv_facets."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from mycnv_facets.annotator import annotate_vcf
from mycnv_facets.cnv_facets import run_cnv_facets
from mycnv_facets.converter import convert_vcf_to_table
from mycnv_facets.normalizer import generate_normal_chrom


@click.command()
# --- Required ---
@click.option("-o", "--out", required=True, help="Output prefix for all output files.")
# --- Input ---
@click.option("-t", "--snp-tumour", default=None, help="BAM file for tumour sample.")
@click.option("-n", "--snp-normal", default=None, help="BAM file for normal sample.")
@click.option("-vcf", "--snp-vcf", default=None, help="VCF file of SNPs for pileup.")
@click.option("-p", "--pileup", default=None, help="Pre-computed pileup CSV file.")
# --- Quality control ---
@click.option("-mq", "--snp-mapq", default=5, type=int, show_default=True, help="Minimum mapping quality threshold.")
@click.option("-bq", "--snp-baq", default=10, type=int, show_default=True, help="Minimum base quality threshold.")
@click.option("-A", "--snp-count-orphans", is_flag=True, help="Do not discard anomalous read pairs.")
@click.option("-N", "--snp-nprocs", default=1, type=int, show_default=True, help="Parallel processes for pileup.")
@click.option("-d", "--depth", default=(25, 4000), type=(int, int), show_default=True, help="Min and max depth in normal sample.")
# --- Analysis ---
@click.option("-T", "--targets", default=None, help="BED file of target regions.")
@click.option("-cv", "--cval", default=(25, 150), type=(float, float), show_default=True, help="Critical values for segmentation.")
@click.option("-snp", "--nbhd-snp", default="auto", show_default=True, help="SNP neighbourhood size or 'auto'.")
@click.option("-a", "--annotation", default=None, help="BED annotation file.")
@click.option("-g", "--gbuild", default="hg19", show_default=True, type=click.Choice(["hg19", "hg38", "mm9", "mm10"]), help="Genome build.")
@click.option("-u", "--unmatched", is_flag=True, help="Treat normal as unmatched.")
@click.option("-np", "--no-cov-plot", is_flag=True, help="Skip coverage plot generation.")
@click.option("-s", "--rnd-seed", default=None, help="Random seed.")
# --- mycnv_facets extensions ---
@click.option("--gene-list", default=None, help="Gene list file (one symbol per line).")
@click.option("--refseq", default=None, help="RefSeq annotation file (ncbiRefSeqCurated format).")
@click.option("--skip-cnv-facets", is_flag=True, help="Skip cnv_facets.R step; start from existing VCF.")
@click.option("--input-vcf", default=None, help="Input CNV VCF file (for --skip-cnv-facets mode).")
@click.option("--skip-annotation", is_flag=True, help="Skip gene annotation step.")
@click.option("--skip-table", is_flag=True, help="Skip VCF-to-table conversion.")
@click.option("--skip-normalizer", is_flag=True, help="Skip normal chromosome generation.")
def main(
    out: str,
    snp_tumour: str | None,
    snp_normal: str | None,
    snp_vcf: str | None,
    pileup: str | None,
    snp_mapq: int,
    snp_baq: int,
    snp_count_orphans: bool,
    snp_nprocs: int,
    depth: tuple[int, int],
    targets: str | None,
    cval: tuple[float, float],
    nbhd_snp: str,
    annotation: str | None,
    gbuild: str,
    unmatched: bool,
    no_cov_plot: bool,
    rnd_seed: str | None,
    gene_list: str | None,
    refseq: str | None,
    skip_cnv_facets: bool,
    input_vcf: str | None,
    skip_annotation: bool,
    skip_table: bool,
    skip_normalizer: bool,
):
    """mycnv_facets — CNV detection, gene annotation, and tabular output.

    Wraps cnv_facets.R (FACETS-based CNV caller) and adds:
    1. Gene panel annotation of CNV intervals
    2. Tabular output (TSV) with simplified fields
    3. Normal chromosome identification
    """
    # Validate inputs
    if not skip_cnv_facets and not pileup:
        if not snp_tumour or not snp_normal:
            click.echo("Error: --snp-tumour and --snp-normal are required unless --skip-cnv-facets or --pileup is set.", err=True)
            sys.exit(1)
        if not snp_vcf:
            click.echo("Error: --snp-vcf is required unless --skip-cnv-facets or --pileup is set.", err=True)
            sys.exit(1)

    # Ensure output directory exists
    out_dir = str(Path(out).parent)
    if out_dir:
        Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Step 1: Run cnv_facets.R
    if not skip_cnv_facets:
        click.echo("[1/4] Running cnv_facets.R...")
        vcf_path = run_cnv_facets(
            tumour_bam=snp_tumour,
            normal_bam=snp_normal,
            snp_vcf=snp_vcf,
            out=out,
            pileup=pileup,
            mapq=snp_mapq,
            baq=snp_baq,
            count_orphans=snp_count_orphans,
            nprocs=snp_nprocs,
            depth=depth,
            targets=targets,
            cval=cval,
            nbhd_snp=nbhd_snp,
            annotation=annotation,
            gbuild=gbuild,
            unmatched=unmatched,
            no_cov_plot=no_cov_plot,
            rnd_seed=rnd_seed,
        )
        click.echo(f"  VCF: {vcf_path}")
    else:
        # --skip-cnv-facets: use the provided VCF directly
        if input_vcf and Path(input_vcf).exists():
            vcf_path = input_vcf
        elif Path(f"{out}.vcf.gz").exists():
            vcf_path = f"{out}.vcf.gz"
        else:
            click.echo(
                f"Error: VCF file not found. Provide --input-vcf or ensure {out}.vcf.gz exists.",
                err=True,
            )
            sys.exit(1)
        click.echo(f"[1/4] Skipping cnv_facets.R, using {vcf_path}")

    # Step 2: Gene annotation
    annotated_vcf = vcf_path
    if gene_list and refseq and not skip_annotation:
        click.echo("[2/4] Annotating VCF with gene panel...")
        annotated_vcf = f"{out}.annotated.vcf.gz"
        annotate_vcf(
            vcf_path=vcf_path,
            out_vcf_path=annotated_vcf,
            gene_list_path=gene_list,
            refseq_path=refseq,
        )
        click.echo(f"  Annotated VCF: {annotated_vcf}")
    elif gene_list or refseq:
        click.echo("Warning: both --gene-list and --refseq required for annotation. Skipping.", err=True)
    else:
        click.echo("[2/4] No gene annotation requested.")

    # Step 3: VCF to TSV
    tsv_path = f"{out}.annotated.tsv"
    if not skip_table:
        click.echo("[3/4] Converting VCF to table...")
        convert_vcf_to_table(
            vcf_path=annotated_vcf,
            out_tsv_path=tsv_path,
        )
        click.echo(f"  TSV: {tsv_path}")
    else:
        click.echo("[3/4] Skipping table conversion.")

    # Step 4: Normal chromosome identification
    if Path(tsv_path).exists() and not skip_normalizer:
        click.echo("[4/4] Generating normal_chrom.txt...")
        result = generate_normal_chrom(
            tsv_path=tsv_path,
            out_path=f"{out}.normal_chrom.txt",
        )
        click.echo(f"  Normal chroms: {result or '(none)'}")
    else:
        click.echo("[4/4] Skipping normal chromosome identification.")

    click.echo("Done.")
