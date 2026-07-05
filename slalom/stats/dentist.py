"""DENTIST-S: a simplified, single-lead-variant form of the DENTIST outlier test.

For a chosen lead variant, DENTIST-S compares each variant's observed Z score against
the value predicted from the lead Z and the variant-to-lead LD (r). A large deviation,
scored as a 1-df chi-square, flags an association-statistic outlier that is inconsistent
with the local LD structure (see Kanai et al. 2022; Chen et al. 2021 for DENTIST).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Tuple

import numpy as np
import scipy.stats

if TYPE_CHECKING:
    import numpy.typing as npt


def dentist_s(z: "npt.ArrayLike", r: "npt.ArrayLike", lead_idx: int) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the DENTIST-S statistic and its -log10 p-value for every variant.

    Parameters
    ----------
    z : array_like
        Marginal Z scores (beta / se) for all variants in the locus.
    r : array_like
        Signed LD (r) between each variant and the lead variant.
    lead_idx : int
        Row index of the lead variant within `z` / `r`.

    Returns
    -------
    t : ndarray
        DENTIST-S test statistic per variant (NaN at the lead variant).
    nlog10p : ndarray
        -log10 p-value under a 1-df chi-square (NaN at the lead variant).
    """
    z = np.asarray(z, dtype=np.float64)
    r = np.asarray(r, dtype=np.float64)
    lead_z = z[lead_idx]

    # r == 1 at the lead gives 0/0 -> NaN (overwritten below); r slightly > 1 from
    # numerical noise gives a negative denominator. Both are handled after the divide,
    # so the expected warnings are silenced here.
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (z - r * lead_z) ** 2 / (1 - r**2)
    # A negative denominator (|r| > 1) makes the statistic meaningless; treat as
    # maximally outlying.
    t = np.where(t < 0, np.inf, t)
    t[lead_idx] = np.nan

    nlog10p = scipy.stats.chi2.logsf(t, df=1) / -np.log(10)
    return t, nlog10p
