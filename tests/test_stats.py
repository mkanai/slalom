import numpy as np

from slalom.stats.abf import abf, get_cs
from slalom.stats.dentist import dentist_s


def test_abf_prob_sums_to_one():
    beta = np.array([0.1, 0.3, -0.2, 0.05])
    se = np.array([0.05, 0.05, 0.05, 0.05])
    lbf, prob = abf(beta, se)
    assert np.isclose(prob.sum(), 1.0)
    # The strongest signal (largest |z|) gets the highest posterior.
    assert prob.argmax() == np.abs(beta / se).argmax()
    assert lbf.shape == beta.shape


def test_get_cs_covers_requested_mass():
    variants = np.array(["v0", "v1", "v2", "v3"])
    prob = np.array([0.6, 0.3, 0.07, 0.03])
    cs = get_cs(variants, prob, coverage=0.95)
    # 0.6 + 0.3 + 0.07 = 0.97 > 0.95, so the top three are needed.
    assert set(cs) == {"v0", "v1", "v2"}


def test_dentist_s_lead_is_nan_and_zero_elsewhere_when_consistent():
    # If every variant's z equals r * lead_z exactly, the statistic is 0 (no outlier).
    lead_z = 5.0
    r = np.array([1.0, 0.5, -0.3, 0.8])
    z = r * lead_z
    t, nlog10p = dentist_s(z, r, lead_idx=0)
    assert np.isnan(t[0]) and np.isnan(nlog10p[0])
    assert np.allclose(t[1:], 0.0)


def test_dentist_s_flags_outlier():
    r = np.array([1.0, 0.0])
    z = np.array([5.0, 5.0])  # variant 1 has high z despite zero LD with the lead
    t, nlog10p = dentist_s(z, r, lead_idx=0)
    assert t[1] == 25.0  # (5 - 0*5)^2 / (1 - 0)
    assert nlog10p[1] > 4
