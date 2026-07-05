"""Lead-variant LD annotation read from a Hail ``BlockMatrix`` via ldcov (no Hail/Spark).

For each locus SLALOM only needs one row of the LD matrix: the correlation (r) between the
lead variant and every other variant. We therefore read just the blocks along the lead
variant's block row/column instead of materialising the full submatrix.

gnomAD LD ``BlockMatrix`` stores are upper-triangular, so r(lead, j) is the stored value at
(min(lead, j), max(lead, j)); the diagonal (j == lead) is 1. Variant-to-matrix-index mapping
is delegated to ldcov's ``VariantIndex.match_variants`` (one region query per contig). By
default a variant must match the panel's exact ref/alt orientation; ``allow_allele_swap=True``
also matches the ref/alt swap and sign-flips r accordingly.
"""

from __future__ import annotations

from collections import defaultdict
from typing import List, Optional

import numpy as np
import pandas as pd
from ldcov.io.blockmatrix import HailBlockMatrixReader
from ldcov.io.variant_index import VariantIndex


def _extract_lead_row(reader: HailBlockMatrixReader, lead_idx: int, target_idxs: List[int]) -> np.ndarray:
    """Return the upper-triangular value r(lead, j) for each j in `target_idxs`.

    Reads only the blocks touched by the lead variant's row/column; off-band pairs that
    are not stored come back as NaN.
    """
    bs = reader.block_size
    out = np.full(len(target_idxs), np.nan, dtype=np.float64)

    # Group requested cells by the block that stores them, so each block is read once.
    groups = defaultdict(list)  # (block_i, block_j) -> [(out_pos, local_i, local_j), ...]
    for pos, g in enumerate(target_idxs):
        a, b = (lead_idx, g) if lead_idx <= g else (g, lead_idx)
        groups[(a // bs, b // bs)].append((pos, a % bs, b % bs))

    for (bi, bj), cells in groups.items():
        block = reader.read_block(bi, bj)
        if block is None:
            continue
        for pos, li, lj in cells:
            out[pos] = block[li, lj]
    return out


def lead_variant_r(
    df: pd.DataFrame,
    lead_row: int,
    bm_path: str,
    variant_index_path: str,
    storage_options: Optional[dict] = None,
    block_cache: int = 8,
    allow_allele_swap: bool = False,
) -> np.ndarray:
    """LD (r) between the lead variant and every row of `df` from one LD panel.

    Parameters
    ----------
    df : DataFrame
        Locus variants with aligned `ref`/`alt` columns (and `chromosome`, `position`).
    lead_row : int
        Row index of the lead variant within `df`.
    bm_path : str
        Hail BlockMatrix directory (local, gs://, or s3://).
    variant_index_path : str
        ldcov Parquet variant index for that matrix.
    allow_allele_swap : bool
        If True, match variants stored in the panel with ref/alt swapped and sign-flip their
        r. Default False matches the exact ref/alt orientation only.

    Returns
    -------
    r : ndarray
        r aligned with `df` rows; NaN where a variant (or the lead) is absent from the panel.
    """
    r = np.full(len(df), np.nan, dtype=np.float64)

    vi = VariantIndex(variant_index_path, storage_options=storage_options)
    query = pd.DataFrame(
        {
            "contig": df["chromosome"].to_numpy(),
            "position": df["position"].to_numpy(),
            "ref": df["ref"].to_numpy(),
            "alt": df["alt"].to_numpy(),
        }
    )
    matches = vi.match_variants(query, allow_allele_swap=allow_allele_swap)
    idx = matches["idx"].to_numpy()
    flip = matches["flip"].to_numpy()

    lead_gi = int(idx[lead_row])
    if lead_gi < 0:
        return r  # lead not in this LD panel: whole column is undefined

    matched = np.where(idx >= 0)[0]
    reader = HailBlockMatrixReader(bm_path, storage_options=storage_options, block_cache=block_cache)
    values = _extract_lead_row(reader, lead_gi, idx[matched].tolist())

    # r changes sign when exactly one of {lead, variant} is stored allele-swapped.
    lead_sign = -1.0 if flip[lead_row] else 1.0
    signs = np.where(flip[matched], -1.0, 1.0) * lead_sign
    r[matched] = values * signs
    return r
