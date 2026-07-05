"""Reference-data constants and default resource paths for SLALOM.

Runtime reads three kinds of reference data, all Hail-free:

1. gnomAD LD ``BlockMatrix`` stores on public GCS (read via ldcov's pure-Python reader).
2. gnomAD LD *variant indices* as Parquet, hosted by ldcov at ``gs://ldcov-requester-pays``
   (map each variant to its BlockMatrix row/column; built once from the ``.ht`` companions).
3. gnomAD sites annotations + conversion-unstable positions (CUPs) as Parquet, built once
   from the source Hail Tables by the helpers under ``scripts/`` (see README).

Every default path here is overridable on the command line, so users can point at their
own converted copies (e.g. in their own requester-pays bucket).
"""

from __future__ import annotations

from typing import List, Optional

# Latest gnomAD release used for annotations, per reference genome build.
GNOMAD_LATEST_VERSIONS = {"GRCh37": "2.1.1", "GRCh38": "3.1.2"}

# Populations with a gnomAD frequency column, per build.
GNOMAD_POPS = {
    "GRCh37": ["afr", "amr", "eas", "fin", "nfe"],
    "GRCh38": ["afr", "amr", "eas", "fin", "nfe", "sas"],
}

# The gnomAD LD BlockMatrices are all GRCh37 (r2.1.1) genomes; only these five
# populations have a released LD matrix. GRCh38 sumstats use the same matrices via a
# lifted-over (b38) variant index.
LD_POPS = ["afr", "amr", "eas", "fin", "nfe"]

# Short build tag used in ldcov's hosted variant-index filenames.
BUILD_TAG = {"GRCh37": "b37", "GRCh38": "b38"}

# Public gnomAD LD BlockMatrix (r values), one per population. Requester-pays-free
# (public bucket), but reads still need GCS credentials configured via gcsfs.
GNOMAD_LD_BM = "gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.{pop}.common.ld.bm"

# ldcov-hosted Parquet variant indices for the matrices above (requester-pays).
# Built from gnomad.genomes.r2.1.1.{pop}.common.adj.ld.variant_indices[.b38].ht.
LDCOV_VARIANT_INDEX_NAME = "gnomad_v2.{pop}.{build}.variant_index.parquet"
LDCOV_VARIANT_INDEX = "gs://ldcov-requester-pays/" + LDCOV_VARIANT_INDEX_NAME

# gnomAD sites annotations (most-severe consequence + per-pop AF) as Parquet, built by
# scripts/make_gnomad_sites_parquet.py. {version} is the build-specific gnomAD version.
GNOMAD_SITES_PARQUET = (
    "gs://finucane-requester-pays/slalom/parquet/" "gnomad.genomes.r{version}.sites.most_severe.{build}.parquet"
)

# Novel conversion-unstable positions (CUPs) + reject regions as Parquet, built by
# scripts/make_cup_parquet.py.
CUP_PARQUET = "gs://finucane-requester-pays/slalom/parquet/FASTA_BED.ALL_{reference_genome}.cups.parquet"

# Consequence categories treated as non-synonymous in the summary.
NONSYN_CONSEQUENCES = ["pLoF", "Missense"]


def gnomad_version(reference_genome: str) -> str:
    return GNOMAD_LATEST_VERSIONS[reference_genome]


def ld_bm_paths() -> List[str]:
    """Default gnomAD LD BlockMatrix paths, one per LD population."""
    return [GNOMAD_LD_BM.format(pop=pop) for pop in LD_POPS]


def ld_variant_index_paths(reference_genome: str, base: Optional[str] = None) -> List[str]:
    """Parquet variant-index paths, one per LD population.

    With ``base=None`` these are the ldcov-hosted requester-pays URLs. Pass a local
    directory (e.g. an extracted ldcov bundle) as ``base`` to resolve local copies for
    offline runs.
    """
    build = BUILD_TAG[reference_genome]
    names = [LDCOV_VARIANT_INDEX_NAME.format(pop=pop, build=build) for pop in LD_POPS]
    if base is None:
        return [LDCOV_VARIANT_INDEX.format(pop=pop, build=build) for pop in LD_POPS]
    base = base.rstrip("/")
    return [f"{base}/{name}" for name in names]


def gnomad_sites_parquet(reference_genome: str) -> str:
    return GNOMAD_SITES_PARQUET.format(version=gnomad_version(reference_genome), build=BUILD_TAG[reference_genome])


def cup_parquet(reference_genome: str) -> str:
    return CUP_PARQUET.format(reference_genome=reference_genome)
