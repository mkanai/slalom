"""gnomAD sites annotations (consequence + frequency) and CUP flags, Hail-free.

gnomAD sites annotations (consequence + frequency) and CUP flags come from region-filtered
Parquet queries plus pandas merges. Join keys use the *aligned* alleles (`ref`/`alt`).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List

import numpy as np
import pandas as pd

from ..resources import GNOMAD_POPS, gnomad_version

if TYPE_CHECKING:
    from ..io.reference import ReferencePanel

_CONSEQUENCE_COLUMNS = ["most_severe", "gene_most_severe", "consequence"]


def _gather_sites(panel: "ReferencePanel", df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """Concatenate per-contig gnomAD sites region queries covering all variants in `df`."""
    frames = []
    for chrom, sub in df.groupby("chromosome"):
        pos = sub["position"]
        frames.append(panel.query_sites(chrom, int(pos.min()), int(pos.max()), columns=columns))
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True)


def annotate_consequence_and_freq(
    df: pd.DataFrame,
    panel: "ReferencePanel",
    reference_genome: str,
    annotate_consequence: bool = False,
    annotate_freq: bool = False,
) -> pd.DataFrame:
    """Left-join gnomAD most-severe consequence and/or per-population AF onto `df`.

    Frequency columns are named ``gnomad_v{major}_af_{pop}`` (e.g. ``gnomad_v3_af_nfe`` for
    GRCh38). Variants absent from gnomAD get NaN/NA.
    """
    if not (annotate_consequence or annotate_freq):
        return df

    version = gnomad_version(reference_genome)
    major = version[0]
    pops = GNOMAD_POPS[reference_genome]

    consequence_cols = _CONSEQUENCE_COLUMNS if annotate_consequence else []
    freq_cols = [f"af_{pop}" for pop in pops] if annotate_freq else []
    columns = ["contig", "position", "ref", "alt"] + consequence_cols + freq_cols

    sites = _gather_sites(panel, df, columns)
    rename = {f"af_{pop}": f"gnomad_v{major}_af_{pop}" for pop in pops}
    sites = sites.rename(columns={"contig": "chromosome", **rename})
    # Guard against duplicate (contig,pos,ref,alt) rows blowing up the row count on merge.
    sites = sites.drop_duplicates(subset=["chromosome", "position", "ref", "alt"])

    df = df.merge(sites, on=["chromosome", "position", "ref", "alt"], how="left")
    return df


def annotate_cups(df: pd.DataFrame, panel: "ReferencePanel") -> pd.DataFrame:
    """Flag variants that fall in a novel CUP or reject interval.

    Adds a boolean `in_cups` column. CUPs are half-open intervals [start, end); a variant
    is flagged when start <= position < end for any interval.
    """
    df = df.copy()
    in_cups = np.zeros(len(df), dtype=bool)
    positions = df["position"].to_numpy()
    chrom = df["chromosome"].astype(str).to_numpy()

    for c in np.unique(chrom):
        rows_i = np.where(chrom == c)[0]
        pos = positions[rows_i]
        cups = panel.query_cups(c, int(pos.min()), int(pos.max()))
        starts = cups["start"].to_numpy()
        ends = cups["end"].to_numpy()
        flags = np.zeros(len(rows_i), dtype=bool)
        for s, e in zip(starts, ends):
            flags |= (pos >= s) & (pos < e)
        in_cups[rows_i] = flags

    df["in_cups"] = in_cups
    return df
