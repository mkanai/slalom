"""Approximate Bayes Factor (ABF) fine-mapping.

Wakefield's ABF gives a per-variant posterior inclusion probability (PIP) and a
credible set from marginal effect sizes and standard errors alone.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Tuple

import numpy as np
import scipy.special

if TYPE_CHECKING:
    import numpy.typing as npt


def abf(beta: "npt.ArrayLike", se: "npt.ArrayLike", W: float = 0.04) -> Tuple[np.ndarray, np.ndarray]:
    """Wakefield approximate Bayes factors and posterior inclusion probabilities.

    Parameters
    ----------
    beta, se : array_like
        Marginal effect size estimates and their standard errors.
    W : float
        Prior variance of the effect size (default 0.04).

    Returns
    -------
    lbf : ndarray
        Log Bayes factor per variant.
    prob : ndarray
        Posterior inclusion probability per variant (sums to 1).
    """
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    z = beta / se
    V = se**2
    r = W / (W + V)
    lbf = 0.5 * (np.log(1 - r) + (r * z**2))
    denom = scipy.special.logsumexp(lbf)
    prob = np.exp(lbf - denom)
    return lbf, prob


def get_cs(variant: "npt.ArrayLike", prob: "npt.ArrayLike", coverage: float = 0.95) -> np.ndarray:
    """Return the credible set of variants covering `coverage` of the posterior mass.

    Returns an empty array if the posterior never reaches `coverage` (e.g. all-NaN PIPs
    from non-finite `se`/`beta`), rather than raising.
    """
    variant = np.asarray(variant)
    prob = np.asarray(prob, dtype=np.float64)
    ordering = np.argsort(prob)[::-1]
    hits = np.where(np.cumsum(prob[ordering]) > coverage)[0]
    if hits.size == 0:
        return variant[:0]
    cs = variant[ordering][: (hits[0] + 1)]
    return cs
