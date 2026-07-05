"""Read the per-locus SNP input and write SLALOM output tables (local or gs://)."""

from __future__ import annotations

from typing import Optional

import fsspec
import pandas as pd

# Columns the SNP file must provide for the core pipeline. Extra columns (e.g. maf,
# per-population sample counts, gamma) are carried through untouched.
REQUIRED_COLUMNS = ["chromosome", "position", "allele1", "allele2", "beta", "se"]


def read_snp(path: str, storage_options: Optional[dict] = None) -> pd.DataFrame:
    """Read a whitespace-delimited SNP file into a DataFrame.

    `chromosome` is read as a string so contig names like "chr1"/"X" survive round-trip.
    """
    with fsspec.open(path, "rt", **(storage_options or {})) as fh:
        df = pd.read_csv(fh, sep=r"\s+", dtype={"chromosome": "string"}, engine="python")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"SNP file {path} is missing required columns: {missing}")
    # Normalise allele columns to plain strings for downstream comparisons.
    df["allele1"] = df["allele1"].astype(str)
    df["allele2"] = df["allele2"].astype(str)
    df["chromosome"] = df["chromosome"].astype(str)
    return df


def write_table(df: pd.DataFrame, path: str, storage_options: Optional[dict] = None) -> None:
    """Write a DataFrame as a tab-separated table with NA-encoded missing values."""
    with fsspec.open(path, "wt", **(storage_options or {})) as fh:
        df.to_csv(fh, sep="\t", na_rep="NA", index=False)
