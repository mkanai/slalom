import numpy as np
import pandas as pd

from slalom.ld import lead_variant_r
from tests.helpers import make_symmetric_bm, write_variant_index


def _setup(tmp_path, n=6):
    full = make_symmetric_bm(str(tmp_path / "ld.bm"), n=n, block_size=3, seed=1)
    positions = [100 * (i + 1) for i in range(n)]
    write_variant_index(
        str(tmp_path / "vi.parquet"),
        contig="1",
        positions=positions,
        refs=["A"] * n,
        alts=["G"] * n,
        idxs=list(range(n)),
    )
    return full, positions


def _df(positions, ref=None, alt=None):
    n = len(positions)
    return pd.DataFrame(
        {
            "chromosome": ["1"] * n,
            "position": positions,
            "ref": ref if ref is not None else ["A"] * n,
            "alt": alt if alt is not None else ["G"] * n,
        }
    )


def test_lead_row_matches_symmetric_matrix(tmp_path):
    full, positions = _setup(tmp_path)
    lead = 2
    r = lead_variant_r(
        _df(positions),
        lead_row=lead,
        bm_path=str(tmp_path / "ld.bm"),
        variant_index_path=str(tmp_path / "vi.parquet"),
    )
    # The extracted row must equal the (symmetrised) lead row of the full matrix.
    assert np.allclose(r, full[lead])
    assert np.isclose(r[lead], 1.0)  # self-correlation on the diagonal


def test_swapped_orientation_is_unmatched_by_default(tmp_path):
    # Default (allow_allele_swap=False): a ref/alt-swapped variant is matched in exact
    # orientation only, so it is NOT found in the panel -> NaN.
    full, positions = _setup(tmp_path)
    ref = ["A"] * len(positions)
    alt = ["G"] * len(positions)
    ref[4], alt[4] = "G", "A"  # swapped relative to the index
    r = lead_variant_r(
        _df(positions, ref, alt),
        lead_row=2,
        bm_path=str(tmp_path / "ld.bm"),
        variant_index_path=str(tmp_path / "vi.parquet"),
    )
    assert np.isnan(r[4])  # swapped variant is unmatched
    assert np.isclose(r[3], full[2, 3])  # exact-orientation entries are unchanged


def test_allow_allele_swap_matches_and_flips_sign(tmp_path):
    # With allow_allele_swap=True the swapped variant matches, and its r is sign-flipped.
    full, positions = _setup(tmp_path)
    ref = ["A"] * len(positions)
    alt = ["G"] * len(positions)
    ref[4], alt[4] = "G", "A"  # swapped relative to the index
    r = lead_variant_r(
        _df(positions, ref, alt),
        lead_row=2,
        bm_path=str(tmp_path / "ld.bm"),
        variant_index_path=str(tmp_path / "vi.parquet"),
        allow_allele_swap=True,
    )
    assert np.isclose(r[4], -full[2, 4])  # matched with a sign flip
    assert np.isclose(r[3], full[2, 3])  # exact-orientation entries unchanged


def test_missing_lead_returns_all_nan(tmp_path):
    _setup(tmp_path)
    # Lead variant absent from the panel (position not in the index) -> whole column NaN.
    df = _df([100], ref=["C"], alt=["T"])  # (100, C, T) not in the (A, G) index
    r = lead_variant_r(
        df,
        lead_row=0,
        bm_path=str(tmp_path / "ld.bm"),
        variant_index_path=str(tmp_path / "vi.parquet"),
    )
    assert np.all(np.isnan(r))
