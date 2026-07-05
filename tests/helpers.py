"""Shared test helpers: synthesize on-disk Hail BlockMatrices and Parquet references.

The BlockMatrix writer mirrors the Hail on-disk format (the same approach ldcov uses in
its own tests) so the LD path can be exercised end-to-end without Hail.
"""

import json
import os
import struct

import lz4.block
import numpy as np
import pandas as pd


def _bm_frame(body: bytes) -> bytes:
    comp = lz4.block.compress(body, store_size=False)
    return struct.pack("<i", len(comp) + 4) + struct.pack("<i", len(body)) + comp


def _bm_encode_block(mat: np.ndarray) -> bytes:
    rows, cols = mat.shape
    body = (
        struct.pack("<i", rows)
        + struct.pack("<i", cols)
        + struct.pack("<b", 0)
        + np.ascontiguousarray(mat, dtype="<f8").tobytes(order="F")
    )
    return _bm_frame(body)


def write_synthetic_bm(path, blocks, n_rows, n_cols, block_size):
    """Write a Hail-format upper-triangular BlockMatrix directory from a block dict."""
    n_block_rows = (n_rows + block_size - 1) // block_size
    items = sorted(blocks.keys(), key=lambda ij: ij[0] + ij[1] * n_block_rows)
    os.makedirs(os.path.join(path, "parts"), exist_ok=True)
    part_files = []
    for slot, ij in enumerate(items):
        pf = "part-%05d" % slot
        part_files.append(pf)
        with open(os.path.join(path, "parts", pf), "wb") as fh:
            fh.write(_bm_encode_block(blocks[ij]))
    meta = {
        "blockSize": block_size,
        "nRows": n_rows,
        "nCols": n_cols,
        "maybeFiltered": [i + j * n_block_rows for (i, j) in items],
        "partFiles": part_files,
    }
    with open(os.path.join(path, "metadata.json"), "w") as fh:
        json.dump(meta, fh)
    return path


def make_symmetric_bm(path, n=6, block_size=3, seed=0):
    """Write an n x n symmetric BM (upper-tri stored). Returns the dense n x n ndarray."""
    rng = np.random.default_rng(seed)
    full = rng.standard_normal((n, n))
    full = (full + full.T) / 2.0
    np.fill_diagonal(full, 1.0)
    nb = (n + block_size - 1) // block_size
    blocks = {}
    for bi in range(nb):
        for bj in range(bi, nb):
            sub = full[bi * block_size : (bi + 1) * block_size, bj * block_size : (bj + 1) * block_size]
            blocks[(bi, bj)] = np.triu(sub).copy() if bi == bj else sub.copy()
    write_synthetic_bm(path, blocks, n, n, block_size)
    return full


def write_variant_index(path, contig, positions, refs, alts, idxs):
    """Write an ldcov-style Parquet variant index."""
    df = pd.DataFrame(
        {
            "contig": [str(contig)] * len(positions),
            "position": list(positions),
            "ref": list(refs),
            "alt": list(alts),
            "idx": list(idxs),
        }
    )
    df.to_parquet(path, index=False)
    return path
