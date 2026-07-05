import pandas as pd

from slalom.annotate.alleles import align_alleles


class FakePanel:
    """Minimal ReferencePanel stand-in: gnomAD variants as a fixed (pos, ref, alt) table."""

    def __init__(self, records):
        self._df = pd.DataFrame(records, columns=["position", "ref", "alt"])

    def query_sites(self, chrom, start, end, columns=None):
        df = self._df[(self._df["position"] >= start) & (self._df["position"] <= end)]
        return df.reset_index(drop=True)


def _df(rows):
    return pd.DataFrame(rows)


def test_align_keeps_matching_orientation():
    panel = FakePanel([(100, "A", "G")])
    df = _df([{"chromosome": "1", "position": 100, "allele1": "A", "allele2": "G", "beta": 0.5}])
    out = align_alleles(df, panel)
    assert out.loc[0, "ref"] == "A" and out.loc[0, "alt"] == "G"
    assert out.loc[0, "beta"] == 0.5  # no flip


def test_align_swaps_and_flips_beta():
    # gnomAD stores G/A; input is A/G -> swap, so beta must be negated.
    panel = FakePanel([(100, "G", "A")])
    df = _df([{"chromosome": "1", "position": 100, "allele1": "A", "allele2": "G", "beta": 0.5}])
    out = align_alleles(df, panel)
    assert out.loc[0, "ref"] == "G" and out.loc[0, "alt"] == "A"
    assert out.loc[0, "beta"] == -0.5


def test_align_strand_flip_no_beta_flip():
    # gnomAD stores T/C (strand flip of A/G); orientation preserved -> beta unchanged.
    panel = FakePanel([(100, "T", "C")])
    df = _df([{"chromosome": "1", "position": 100, "allele1": "A", "allele2": "G", "beta": 0.5}])
    out = align_alleles(df, panel)
    assert out.loc[0, "ref"] == "T" and out.loc[0, "alt"] == "C"
    assert out.loc[0, "beta"] == 0.5


def test_align_no_match_keeps_observed():
    panel = FakePanel([(100, "A", "G")])
    df = _df([{"chromosome": "1", "position": 200, "allele1": "C", "allele2": "T", "beta": 0.5}])
    out = align_alleles(df, panel)
    assert out.loc[0, "ref"] == "C" and out.loc[0, "alt"] == "T"
    assert out.loc[0, "beta"] == 0.5
