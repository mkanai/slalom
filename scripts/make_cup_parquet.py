#!/usr/bin/env python
"""Offline builder: novel-CUP + reject interval .ht -> Parquet for SLALOM's --annotate-cups.

Run ONCE per reference genome on a machine with Hail installed. The SLALOM runtime reads
the resulting Parquet with pyarrow (no Hail). The CUP tables are *interval*-keyed (BED
regions), so the output stores half-open intervals: columns contig, start, end, sorted by
(contig, start). A variant is "in CUPs" when its position falls in any interval, matching
the original ``is_defined(cup[locus]) | is_defined(reject[locus])`` logic.

Usage:
    python scripts/make_cup_parquet.py \
        --cup-ht    gs://.../FASTA_BED.ALL_GRCh38.novel_CUPs.ht \
        --reject-ht gs://.../FASTA_BED.ALL_GRCh38.reject_2.ht \
        --out FASTA_BED.ALL_GRCh38.cups.parquet
"""

import argparse


def _interval_frame(ht):  # pragma: no cover (needs Hail)
    import pandas as pd

    ht = ht.key_by()
    ht = ht.select(
        contig=ht.interval.start.contig,
        start=ht.interval.start.position,
        end=ht.interval.end.position,
    )
    rows = ht.collect()
    return pd.DataFrame(
        {
            "contig": [str(r.contig) for r in rows],
            "start": [int(r.start) for r in rows],
            "end": [int(r.end) for r in rows],
        }
    )


def build(cup_ht_path, reject_ht_path, out_path, compression="zstd"):  # pragma: no cover (needs Hail)
    import hail as hl
    import pandas as pd

    hl.init()
    frames = [_interval_frame(hl.read_table(cup_ht_path))]
    if reject_ht_path:
        frames.append(_interval_frame(hl.read_table(reject_ht_path)))

    df = pd.concat(frames, ignore_index=True)
    df = df.drop_duplicates(subset=["contig", "start", "end"])
    df = df.sort_values(["contig", "start"]).reset_index(drop=True)
    # Enforce dtypes so an empty result still writes a string/int64 schema (an all-empty
    # frame would otherwise infer float64 columns and break contig == "chrN" filters).
    df = df.astype({"contig": "string", "start": "int64", "end": "int64"})
    df.to_parquet(out_path, index=False, compression=compression)
    print(f"Wrote {len(df)} CUP intervals to {out_path}")


def main():  # pragma: no cover
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cup-ht", required=True, help="Novel-CUP interval .ht (local/gs://)")
    ap.add_argument("--reject-ht", help="Reject interval .ht (local/gs://); optional")
    ap.add_argument("--out", required=True, help="Output Parquet path")
    ap.add_argument("--compression", default="zstd", help="Parquet compression codec")
    args = ap.parse_args()
    build(args.cup_ht, args.reject_ht, args.out, compression=args.compression)


if __name__ == "__main__":
    main()
