"""Align SNP alleles to the gnomAD reference orientation (Hail-free).

For each variant we try the observed orientation, the ref/alt swap, the strand flip, and the
strand-flip + swap, and adopt the first orientation that exists in gnomAD. When the adopted
orientation swaps ref/alt, the effect size is negated so `beta` stays on the (aligned) alt
allele.

The aligned alleles are written to `ref`/`alt` (used for LD matching and the gnomAD
annotation join); the input `allele1`/`allele2` columns are left untouched so output variant
IDs match the input.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    import pandas as pd

    from ..io.reference import ReferencePanel

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _flip_strand(allele: str) -> str:
    """Strand-complement a single-nucleotide allele; leave anything else unchanged.

    Only the four bases are complemented; indels and multi-base alleles pass through as-is.
    """
    return _COMPLEMENT.get(allele, allele)


def align_alleles(
    df: "pd.DataFrame", panel: "ReferencePanel", reference_genome: Optional[str] = None
) -> "pd.DataFrame":
    """Return `df` with aligned `ref`/`alt` columns and strand/swap-corrected `beta`.

    Parameters
    ----------
    df : DataFrame
        Must contain chromosome, position, allele1, allele2, beta.
    panel : ReferencePanel
        Provides gnomAD sites for the locus (used as an allele-existence oracle).
    reference_genome : str, optional
        Accepted for signature symmetry with the other annotators; unused.
    """
    df = df.reset_index(drop=True)
    n = len(df)
    chrom = df["chromosome"].astype(str).to_numpy()
    position = df["position"].to_numpy()
    a1 = df["allele1"].astype(str).to_numpy(dtype=object)
    a2 = df["allele2"].astype(str).to_numpy(dtype=object)

    ref = a1.copy()
    alt = a2.copy()
    flip = np.zeros(n, dtype=bool)

    for c in np.unique(chrom):
        rows_i = np.where(chrom == c)[0]
        pos = position[rows_i]
        sites = panel.query_sites(c, int(pos.min()), int(pos.max()), columns=["position", "ref", "alt"])
        present = set(zip(sites["position"].tolist(), sites["ref"].tolist(), sites["alt"].tolist()))
        for i in rows_i:
            p = int(position[i])
            b1, b2 = str(a1[i]), str(a2[i])
            f1, f2 = _flip_strand(b1), _flip_strand(b2)
            # (aligned ref, aligned alt, ref/alt swapped?) tried in priority order.
            for r, a, swapped in ((b1, b2, False), (b2, b1, True), (f1, f2, False), (f2, f1, True)):
                if (p, r, a) in present:
                    ref[i], alt[i], flip[i] = r, a, swapped
                    break
            # default (no match): keep observed orientation, no flip (already set).

    df["ref"] = ref.astype(str)
    df["alt"] = alt.astype(str)
    df.loc[flip, "beta"] = -df.loc[flip, "beta"]
    return df
