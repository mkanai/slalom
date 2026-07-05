import numpy as np
import pandas as pd
import pytest

from slalom.pipeline import SlalomConfig, run_slalom
from tests.helpers import make_symmetric_bm, write_variant_index


def _build_locus(tmp_path):
    n = 6
    full = make_symmetric_bm(str(tmp_path / "ld.bm"), n=n, block_size=3, seed=2)
    positions = [100 * (i + 1) for i in range(n)]
    write_variant_index(str(tmp_path / "vi.parquet"), "1", positions, ["A"] * n, ["G"] * n, list(range(n)))

    # gnomAD sites: all six variants, with one missense to exercise the nonsyn summary.
    sites = pd.DataFrame(
        {
            "contig": ["1"] * n,
            "position": positions,
            "ref": ["A"] * n,
            "alt": ["G"] * n,
            "most_severe": ["missense_variant"] * n,
            "gene_most_severe": ["GENE"] * n,
            "consequence": ["Missense", "synonymous", "synonymous", "synonymous", "synonymous", "synonymous"],
            "af_afr": [0.1] * n,
            "af_amr": [0.1] * n,
            "af_eas": [0.1] * n,
            "af_fin": [0.1] * n,
            "af_nfe": [0.1] * n,
        }
    )
    sites.to_parquet(str(tmp_path / "sites.parquet"), index=False)

    # CUPs are half-open intervals [start, end); [300, 301) covers the variant at pos 300.
    cups = pd.DataFrame({"contig": ["1"], "start": [300], "end": [301]})
    cups.to_parquet(str(tmp_path / "cups.parquet"), index=False)

    snp = pd.DataFrame(
        {
            "rsid": [f"rs{i}" for i in range(n)],
            "chromosome": ["1"] * n,
            "position": positions,
            "allele1": ["A"] * n,
            "allele2": ["G"] * n,
            "beta": [0.2, 0.5, -0.1, 0.3, 0.05, 0.4],
            "se": [0.05] * n,
            "p": [1e-3, 1e-8, 0.2, 1e-2, 0.5, 1e-4],
            "n_samples": [10000] * n,
            "n_cases": [2000] * n,
        }
    )
    snp_path = str(tmp_path / "locus.snp")
    snp.to_csv(snp_path, sep="\t", index=False)
    return snp_path, full, positions


def _config(tmp_path, snp_path, **overrides):
    cfg = dict(
        snp=snp_path,
        out=str(tmp_path / "out.txt"),
        reference_genome="GRCh37",
        ld_reference="custom",
        custom_ld_path=str(tmp_path / "ld.bm"),
        custom_ld_variant_index_path=str(tmp_path / "vi.parquet"),
        custom_ld_label="test",
        export_r=True,
        gnomad_sites_parquet=str(tmp_path / "sites.parquet"),
        cup_parquet=str(tmp_path / "cups.parquet"),
    )
    cfg.update(overrides)
    return SlalomConfig(**cfg)


def test_end_to_end_custom_ld(tmp_path):
    snp_path, full, positions = _build_locus(tmp_path)
    cfg = _config(
        tmp_path,
        snp_path,
        align_alleles=True,
        annotate_cups=True,
        annotate_consequence=True,
        annotate_gnomad_freq=True,
        abf=True,
        dentist_s=True,
        summary=True,
        case_control=True,
        lead_variant_choice="p",
    )
    df = run_slalom(cfg)

    # Lead is the min-p variant (index 1).
    lead_idx = 1
    assert bool(df.loc[lead_idx, "lead_variant"]) is True

    # r column equals the custom-panel lead row; r matches the symmetric matrix.
    assert np.allclose(df["r"].to_numpy(), full[lead_idx])
    assert np.isclose(df.loc[lead_idx, "r"], 1.0)

    # DENTIST-S is NaN at the lead and finite elsewhere.
    assert np.isnan(df.loc[lead_idx, "t_dentist_s"])

    # Annotations landed.
    assert bool(df.loc[2, "in_cups"]) is True  # position 300
    assert "gnomad_v2_af_nfe" in df.columns
    assert df.loc[0, "consequence"] == "Missense"

    # Output file exists and omits internal-only columns.
    out = pd.read_csv(cfg.out, sep="\t")
    for internal in ("variant", "ref", "alt"):
        assert internal not in out.columns
    assert "test_lead_r" in out.columns

    # Summary written with expected fields.
    summary = pd.read_csv(cfg.out_summary, sep="\t")
    assert summary.loc[0, "n_total"] == len(positions)
    assert "n_nonsyn" in summary.columns
    assert "min_neff_r2" in summary.columns


def test_default_ref_alt_without_align(tmp_path):
    snp_path, full, _ = _build_locus(tmp_path)
    cfg = _config(tmp_path, snp_path, dentist_s=True, lead_variant_choice="p")
    df = run_slalom(cfg)
    # Without --align-alleles, ref/alt fall back to the observed alleles and LD still resolves.
    assert np.allclose(df["r"].to_numpy(), full[1])


def test_custom_ld_without_export_r_emits_r2(tmp_path):
    # Default export_r=False must not crash and must keep signed r for DENTIST-S while
    # emitting an r^2 column (regression: r/r2 column mismatch fed r^2 into DENTIST-S).
    snp_path, full, _ = _build_locus(tmp_path)
    cfg = _config(tmp_path, snp_path, export_r=False, dentist_s=True, lead_variant_choice="p")
    df = run_slalom(cfg)
    out = pd.read_csv(cfg.out, sep="\t")
    assert "test_lead_r2" in out.columns  # r^2 output column
    assert "test_lead_r" not in out.columns  # signed-r internal column dropped
    # df["r"] stays signed (combined), so DENTIST-S sees signed r, not r^2.
    assert np.allclose(df["r"].to_numpy(), full[1])
    assert np.allclose(out["test_lead_r2"].to_numpy(), full[1] ** 2)


def test_config_validates_gnomad_defaults():
    # A minimal gnomAD-reference config is valid and defaults out_summary from out.
    cfg = SlalomConfig(snp="in.snp", out="out.txt")
    assert cfg.out_summary == "out.summary.txt"


def test_config_custom_requires_paths():
    # ld_reference="custom" without the three custom_ld_* fields is rejected at construction,
    # so the library API is as safe as the CLI (not only cli.main).
    with pytest.raises(ValueError, match="custom_ld_path"):
        SlalomConfig(snp="in.snp", out="out.txt", ld_reference="custom")


def test_config_weighted_average_requires_export_r():
    with pytest.raises(ValueError, match="export_r"):
        SlalomConfig(snp="in.snp", out="out.txt", weighted_average_r={"nfe": "n_nfe"})


def test_config_summary_requires_dentist_s_and_abf():
    with pytest.raises(ValueError, match="dentist_s"):
        SlalomConfig(snp="in.snp", out="out.txt", summary=True, abf=True)
    with pytest.raises(ValueError, match="abf"):
        SlalomConfig(snp="in.snp", out="out.txt", summary=True, dentist_s=True)


def test_gnomad_default_export_r_false_combines_r(tmp_path, monkeypatch):
    # Regression: the default gnomAD path (export_r=False) must combine df['r'] from the
    # signed per-pop columns, not KeyError on a missing gnomad_lead_r_nfe.
    from slalom import pipeline

    snp_path, full, _ = _build_locus(tmp_path)
    n = len(full)

    # Route every gnomAD population to the one synthetic BM/index so no network is touched.
    monkeypatch.setattr(pipeline.resources, "ld_bm_paths", lambda: [str(tmp_path / "ld.bm")] * 5)
    monkeypatch.setattr(
        pipeline.resources,
        "ld_variant_index_paths",
        lambda rg, base=None: [str(tmp_path / "vi.parquet")] * 5,
    )

    cfg = _config(tmp_path, snp_path, ld_reference="gnomad", export_r=False, dentist_s=True, lead_variant_choice="p")
    df = run_slalom(cfg)
    out = pd.read_csv(cfg.out, sep="\t")
    assert "gnomad_lead_r2_nfe" in out.columns
    assert "gnomad_lead_r_nfe" not in out.columns
    # nfe is the default combine source; df['r'] is signed.
    assert np.allclose(df["r"].to_numpy(), full[1])
    assert not np.isnan(df.loc[0, "t_dentist_s"]) or df.loc[0, "lead_variant"]
    assert n == 6
