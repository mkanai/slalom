"""``slalom`` command-line entry point."""

import argparse
import json
import logging
import os
import sys

from .. import __version__, resources
from ..pipeline import SlalomConfig, run_slalom

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger("slalom")


class ParseKwargs(argparse.Action):
    """Parse ``key=value`` pairs into a dict, coercing numeric values to float."""

    def __call__(self, parser, namespace, values, option_string=None):
        result = {}
        for value in values:
            key, val = value.split("=")
            if val.replace(".", "", 1).isnumeric():
                val = float(val)
            result[key] = val
        setattr(namespace, self.dest, result)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="slalom",
        description="SLALOM: flag suspicious loci for meta-analysis fine-mapping using "
        "association-statistic outliers against a local LD reference.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    io = parser.add_argument_group("input / output")
    io.add_argument("--snp", required=True, help="Input per-locus SNP file from fine-mapping")
    io.add_argument("--out", required=True, help="Output table path (local or gs://)")
    io.add_argument("--out-summary", help="Output summary path (default: <out> with .summary.txt suffix)")
    io.add_argument(
        "--reference-genome",
        default="GRCh37",
        choices=["GRCh37", "GRCh38"],
        help="Reference genome of the summary statistics",
    )
    io.add_argument(
        "--storage-options",
        help="JSON dict of fsspec storage options for remote reads/writes "
        '(e.g. \'{"requester_pays": true, "project": "my-gcp-project"}\')',
    )

    lead = parser.add_argument_group("lead variant")
    lead.add_argument("--lead-variant", help="Lead variant (chrom:pos:ref:alt) to anchor LD")
    lead.add_argument(
        "--lead-variant-choice",
        default="p",
        choices=["p", "prob", "gamma", "gamma-p"],
        help="Strategy for choosing the lead variant when --lead-variant is not given",
    )

    ann = parser.add_argument_group("annotations (Parquet reference tables)")
    ann.add_argument("--align-alleles", action="store_true", help="Align alleles to the gnomAD orientation")
    ann.add_argument(
        "--annotate-cups",
        action="store_true",
        help="Annotate novel conversion-unstable positions (CUPs)",
    )
    ann.add_argument("--annotate-consequence", action="store_true", help="Annotate VEP most-severe consequence")
    ann.add_argument("--annotate-gnomad-freq", action="store_true", help="Annotate gnomAD population frequencies")
    ann.add_argument("--gnomad-sites-parquet", help="Override path to the gnomAD sites Parquet")
    ann.add_argument("--cup-parquet", help="Override path to the CUP Parquet")

    ld = parser.add_argument_group("LD reference")
    ld.add_argument("--ld-reference", default="gnomad", choices=["gnomad", "custom"], help="LD reference source")
    ld.add_argument("--custom-ld-path", help="Path to a user-provided LD BlockMatrix")
    ld.add_argument("--custom-ld-variant-index-path", help="Path to a user-provided Parquet variant index")
    ld.add_argument("--custom-ld-label", help="Column label prefix for the user-provided LD")
    ld.add_argument(
        "--ld-variant-index-dir",
        help="Local directory of gnomAD LD variant-index Parquets (e.g. an extracted ldcov "
        "bundle) to use instead of the hosted defaults, avoiding requester-pays reads",
    )
    ld.add_argument("--export-r", action="store_true", help="Export signed r values instead of r^2")
    ld.add_argument(
        "--weighted-average-r",
        nargs="+",
        action=ParseKwargs,
        metavar="POP=WEIGHT",
        help="Weighted-average r across populations, e.g. nfe=n_nfe afr=n_afr "
        "(weight may be a column name or a constant); requires --export-r",
    )
    ld.add_argument(
        "--block-cache",
        type=int,
        default=8,
        help="Number of decoded BlockMatrix blocks to cache per LD panel",
    )

    stat = parser.add_argument_group("statistics")
    stat.add_argument("--abf", action="store_true", help="Run ABF fine-mapping (PIP + credible set)")
    stat.add_argument("--abf-prior-variance", type=float, default=0.04, help="ABF prior effect-size variance")
    stat.add_argument("--dentist-s", action="store_true", help="Annotate the DENTIST-S statistic")

    summ = parser.add_argument_group("summary")
    summ.add_argument("--summary", action="store_true", help="Write a per-locus summary file")
    summ.add_argument("--case-control", action="store_true", help="Treat the input as a case-control study")
    summ.add_argument("--r2-threshold", type=float, default=0.6, help="r^2 threshold for DENTIST-S outliers")
    summ.add_argument(
        "--nlog10p-dentist-s-threshold",
        type=float,
        default=4.0,
        help="-log10 DENTIST-S p-value threshold for outliers",
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser


def _parse_storage_options(raw):
    if not raw:
        return None
    try:
        opts = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"--storage-options is not valid JSON: {exc}") from exc
    if not isinstance(opts, dict):
        raise ValueError("--storage-options must be a JSON object (dict)")
    return opts


def _validate(args):
    if args.ld_reference == "custom" and (
        not args.custom_ld_path or not args.custom_ld_variant_index_path or not args.custom_ld_label
    ):
        raise ValueError(
            "--ld-reference custom requires --custom-ld-path, " "--custom-ld-variant-index-path, and --custom-ld-label"
        )
    if args.weighted_average_r is not None and not args.export_r:
        raise ValueError("--weighted-average-r requires --export-r (it averages signed r).")
    if args.summary and not args.dentist_s:
        raise ValueError("--summary requires --dentist-s (it reports DENTIST-S outliers).")
    if args.summary and not args.abf:
        raise ValueError("--summary requires --abf (it reports max PIP / credible sets).")


def _config_from_args(args):
    ld_variant_index_paths = None
    if args.ld_variant_index_dir:
        ld_variant_index_paths = resources.ld_variant_index_paths(args.reference_genome, base=args.ld_variant_index_dir)

    return SlalomConfig(
        snp=args.snp,
        out=args.out,
        out_summary=args.out_summary,
        reference_genome=args.reference_genome,
        lead_variant=args.lead_variant,
        lead_variant_choice=args.lead_variant_choice,
        align_alleles=args.align_alleles,
        annotate_cups=args.annotate_cups,
        annotate_consequence=args.annotate_consequence,
        annotate_gnomad_freq=args.annotate_gnomad_freq,
        ld_reference=args.ld_reference,
        custom_ld_path=args.custom_ld_path,
        custom_ld_variant_index_path=args.custom_ld_variant_index_path,
        custom_ld_label=args.custom_ld_label,
        export_r=args.export_r,
        weighted_average_r=args.weighted_average_r,
        dentist_s=args.dentist_s,
        abf=args.abf,
        abf_prior_variance=args.abf_prior_variance,
        summary=args.summary,
        case_control=args.case_control,
        r2_threshold=args.r2_threshold,
        nlog10p_dentist_s_threshold=args.nlog10p_dentist_s_threshold,
        gnomad_sites_parquet=args.gnomad_sites_parquet,
        cup_parquet=args.cup_parquet,
        ld_variant_index_paths=ld_variant_index_paths,
        storage_options=_parse_storage_options(args.storage_options),
        block_cache=args.block_cache,
    )


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.verbose:
        logging.getLogger("slalom").setLevel(logging.DEBUG)

    try:
        _validate(args)
        cfg = _config_from_args(args)
    except ValueError as exc:
        parser.error(str(exc))

    # Create the local output directory if needed (gs:// paths are created implicitly).
    if not args.out.startswith(("gs://", "s3://")):
        out_dir = os.path.dirname(os.path.abspath(args.out))
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    try:
        run_slalom(cfg)
    except Exception as exc:  # surface a clean message, non-zero exit for pipelines
        logger.error("Error: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
