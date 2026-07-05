#!/usr/bin/env python
"""Offline builder: gnomAD sites most_severe .ht -> Parquet for SLALOM annotations.

Run ONCE per reference genome on a machine with Hail installed. The SLALOM runtime never
imports Hail; it reads the resulting Parquet with pyarrow (see slalom/io/reference.py).

Output columns: contig, position, ref, alt, most_severe, gene_most_severe, consequence,
and one af_<pop> column per population. The genome-wide sites HTs are large (gnomAD v3 has
~760M variants across ~115k partitions), so the table is written via Spark (map-only, no
shuffle) as a *partitioned Parquet directory*, coalesced to a modest file count. The HT is
already sorted by locus, so each output file covers a contiguous locus range and pyarrow
predicate pushdown on (contig, position) stays efficient. pyarrow.dataset reads the
directory transparently, so point --gnomad-sites-parquet at it like any Parquet path.

Usage:
    python scripts/make_gnomad_sites_parquet.py \
        --ht gs://.../gnomad.genomes.r3.1.2.sites.most_severe.ht \
        --reference-genome GRCh38 \
        --out gs://YOUR_BUCKET/.../gnomad.genomes.r3.1.2.sites.most_severe.b38.parquet
"""

import argparse

# Populations mirror slalom.resources.GNOMAD_POPS; kept local so this script has no
# runtime import of the slalom package (which need not be installed on the Hail machine).
GNOMAD_POPS = {
    "GRCh37": ["afr", "amr", "eas", "fin", "nfe"],
    "GRCh38": ["afr", "amr", "eas", "fin", "nfe", "sas"],
}

CONSEQUENCE_COLUMNS = ["most_severe", "gene_most_severe", "consequence"]


def build(ht_path, out_path, reference_genome, coalesce=256):  # pragma: no cover (needs Hail)
    import hail as hl

    hl.init()
    ht = hl.read_table(ht_path)
    pops = GNOMAD_POPS[reference_genome]

    ht = ht.annotate(
        contig=ht.locus.contig,
        position=ht.locus.position,
        ref=ht.alleles[0],
        alt=ht.alleles[1],
    )
    # Build freq_expr against the post-key_by table so all select fields share one source.
    ht = ht.key_by()
    freq_expr = {f"af_{pop}": ht.freq[pop].AF for pop in pops}
    ht = ht.select("contig", "position", "ref", "alt", *CONSEQUENCE_COLUMNS, **freq_expr)

    # Spark write scales to the full HT; coalesce (no shuffle) trims the output file count.
    ht.to_spark().coalesce(coalesce).write.mode("overwrite").parquet(out_path)
    print(f"Wrote sites Parquet ({len(ht.row)} cols, coalesce={coalesce}) to {out_path}")


def main():  # pragma: no cover
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ht", required=True, help="gnomAD sites most_severe .ht (local/gs://)")
    ap.add_argument("--reference-genome", required=True, choices=["GRCh37", "GRCh38"], help="Build of the .ht")
    ap.add_argument("--out", required=True, help="Output Parquet directory (local or gs://)")
    ap.add_argument(
        "--coalesce",
        type=int,
        default=256,
        help="Number of output part files (fewer = larger files; ~64 for GRCh37, ~256 for GRCh38)",
    )
    args = ap.parse_args()
    build(args.ht, args.out, args.reference_genome, coalesce=args.coalesce)


if __name__ == "__main__":
    main()
