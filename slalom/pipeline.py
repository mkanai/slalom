"""End-to-end SLALOM pipeline: annotate a locus, score DENTIST-S outliers, summarise.

LD is read from a Hail ``BlockMatrix`` through ldcov; annotations come from Parquet
reference tables. No Hail/Spark is required at runtime.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union, cast

import numpy as np

from . import resources
from .annotate.alleles import align_alleles
from .annotate.annotations import annotate_consequence_and_freq, annotate_cups
from .io.reference import ReferencePanel
from .io.snp import read_snp, write_table
from .ld import lead_variant_r
from .stats.abf import abf, get_cs
from .stats.dentist import dentist_s

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class SlalomConfig:
    """Configuration for a single-locus SLALOM run (mirrors the CLI flags)."""

    snp: str
    out: str
    out_summary: Optional[str] = None
    reference_genome: str = "GRCh37"

    lead_variant: Optional[str] = None
    lead_variant_choice: str = "p"

    align_alleles: bool = False
    annotate_cups: bool = False
    annotate_consequence: bool = False
    annotate_gnomad_freq: bool = False

    ld_reference: str = "gnomad"
    custom_ld_path: Optional[str] = None
    custom_ld_variant_index_path: Optional[str] = None
    custom_ld_label: Optional[str] = None
    export_r: bool = False
    weighted_average_r: Optional[Dict[str, Union[str, float]]] = None

    dentist_s: bool = False
    abf: bool = False
    abf_prior_variance: float = 0.04

    summary: bool = False
    case_control: bool = False
    r2_threshold: float = 0.6
    nlog10p_dentist_s_threshold: float = 4.0

    # Reference-data locations (default to the hosted copies in resources.py).
    gnomad_sites_parquet: Optional[str] = None
    cup_parquet: Optional[str] = None
    ld_variant_index_paths: Optional[list] = None
    ld_bm_paths: Optional[list] = None

    # fsspec storage options for remote reads/writes (e.g. requester-pays project).
    storage_options: Optional[dict] = None
    block_cache: int = 8

    def __post_init__(self) -> None:
        if self.out_summary is None:
            self.out_summary = f"{os.path.splitext(self.out)[0]}.summary.txt"


def _make_variant_ids(df: "pd.DataFrame") -> "pd.Series":
    """chrom:pos:ref:alt using the aligned alleles (matches hl.variant_str)."""
    return (
        df["chromosome"].astype(str)
        + ":"
        + df["position"].astype(str)
        + ":"
        + df["ref"].astype(str)
        + ":"
        + df["alt"].astype(str)
    )


def _choose_lead_index(df: "pd.DataFrame", cfg: SlalomConfig) -> int:
    """Return the integer row index of the lead variant per the configured strategy."""
    if cfg.lead_variant is not None:
        matches = np.where(df["variant"].to_numpy() == cfg.lead_variant)[0]
        if len(matches) == 0:
            raise ValueError(f"Lead variant {cfg.lead_variant} not found in the input.")
        return int(matches[0])

    if cfg.lead_variant_choice == "p":
        return int(df["p"].idxmin())
    if cfg.lead_variant_choice == "prob":
        return int(df["prob"].idxmax())
    if cfg.lead_variant_choice in ("gamma", "gamma-p"):
        gamma_idx = np.where(df["gamma"].to_numpy())[0]
        if len(gamma_idx) == 0:
            if cfg.lead_variant_choice == "gamma-p":
                return int(df["p"].idxmin())
            raise ValueError("No lead variants found with gamma.")
        if len(gamma_idx) > 1:
            raise ValueError("Multiple lead variants found with gamma.")
        return int(gamma_idx[0])
    raise ValueError(f"Unknown lead-variant choice: {cfg.lead_variant_choice}")


def _ld_targets(cfg: SlalomConfig) -> Tuple[List[str], List[str], List[str], List[str]]:
    """Return (bm_paths, variant_index_paths, r_labels, r2_labels) for the LD reference.

    ``r_labels`` are the canonical signed-r columns consumed internally (``_combine_ld`` and
    DENTIST-S); ``r2_labels`` are the paired r^2 output names emitted when ``export_r`` is off.
    """
    if cfg.ld_reference == "gnomad":
        bm_paths = cfg.ld_bm_paths or resources.ld_bm_paths()
        index_paths = cfg.ld_variant_index_paths or resources.ld_variant_index_paths(cfg.reference_genome)
        r_labels = [f"gnomad_lead_r_{pop}" for pop in resources.LD_POPS]
        r2_labels = [f"gnomad_lead_r2_{pop}" for pop in resources.LD_POPS]
        return bm_paths, index_paths, r_labels, r2_labels
    # custom single-panel reference (paths are validated non-None for --ld-reference custom)
    r_labels = [f"{cfg.custom_ld_label}_lead_r"]
    r2_labels = [f"{cfg.custom_ld_label}_lead_r2"]
    bm_paths = cast(List[str], [cfg.custom_ld_path])
    index_paths = cast(List[str], [cfg.custom_ld_variant_index_path])
    return bm_paths, index_paths, r_labels, r2_labels


def run_slalom(cfg: SlalomConfig) -> "pd.DataFrame":
    """Run the full SLALOM pipeline for one locus and write the output table(s).

    Returns the annotated per-variant DataFrame.
    """
    panel = ReferencePanel(
        sites_path=cfg.gnomad_sites_parquet or resources.gnomad_sites_parquet(cfg.reference_genome),
        cup_path=cfg.cup_parquet or resources.cup_parquet(cfg.reference_genome),
        storage_options=cfg.storage_options,
    )

    df = read_snp(cfg.snp, storage_options=cfg.storage_options)
    df = df.reset_index(drop=True)

    # Aligned alleles (ref/alt) drive LD matching and the gnomAD annotation join. Without
    # --align-alleles they are just the observed alleles.
    if cfg.align_alleles:
        df = align_alleles(df, panel, cfg.reference_genome)
    else:
        df["ref"] = df["allele1"].astype(str)
        df["alt"] = df["allele2"].astype(str)

    if cfg.annotate_cups:
        df = annotate_cups(df, panel)

    if cfg.annotate_consequence or cfg.annotate_gnomad_freq:
        df = annotate_consequence_and_freq(
            df,
            panel,
            cfg.reference_genome,
            annotate_consequence=cfg.annotate_consequence,
            annotate_freq=cfg.annotate_gnomad_freq,
        )

    df["variant"] = _make_variant_ids(df)

    if cfg.abf:
        lbf, prob = abf(df["beta"], df["se"], W=cfg.abf_prior_variance)
        df["lbf"] = lbf
        df["prob"] = prob
        cs = get_cs(df["variant"].to_numpy(), prob, coverage=0.95)
        cs_99 = get_cs(df["variant"].to_numpy(), prob, coverage=0.99)
        df["cs"] = df["variant"].isin(cs)
        df["cs_99"] = df["variant"].isin(cs_99)

    lead_idx = _choose_lead_index(df, cfg)
    cfg.lead_variant = df["variant"].iloc[lead_idx]
    df["lead_variant"] = False
    df.iloc[lead_idx, df.columns.get_loc("lead_variant")] = True

    # Annotate LD between the lead variant and every variant. Columns hold signed r
    # internally (needed by _combine_ld and DENTIST-S); r^2 is derived for output below.
    bm_paths, index_paths, r_labels, r2_labels = _ld_targets(cfg)
    for bm_path, index_path, col in zip(bm_paths, index_paths, r_labels):
        logger.info("Reading LD for %s from %s", col, bm_path)
        df[col] = lead_variant_r(
            df,
            lead_row=lead_idx,
            bm_path=bm_path,
            variant_index_path=index_path,
            storage_options=cfg.storage_options,
            block_cache=cfg.block_cache,
            # --align-alleles is a request to reconcile orientation with the reference, so it
            # also enables ref/alt-swap matching against the LD panel.
            allow_allele_swap=cfg.align_alleles,
        )

    _combine_ld(df, cfg, r_labels)

    if cfg.dentist_s:
        z = (df["beta"] / df["se"]).to_numpy()
        t, nlog10p = dentist_s(z, df["r"].to_numpy(), lead_idx)
        df["t_dentist_s"] = t
        df["nlog10p_dentist_s"] = nlog10p

    # Unless exporting signed r, convert the per-panel LD columns to r^2 for output.
    if not cfg.export_r:
        for r_col, r2_col in zip(r_labels, r2_labels):
            df[r2_col] = df[r_col] ** 2
        df = df.drop(columns=list(r_labels))

    # Output table: drop internal-only columns.
    out_df = df.drop(columns=["variant", "ref", "alt"])
    write_table(out_df, cfg.out, storage_options=cfg.storage_options)
    logger.info("Wrote %s", cfg.out)

    if cfg.summary:
        summary_df = _build_summary(df, cfg)
        # out_summary is always populated in __post_init__ (defaulted from `out`).
        write_table(summary_df, cast(str, cfg.out_summary), storage_options=cfg.storage_options)
        logger.info("Wrote %s", cfg.out_summary)

    return df


def _combine_ld(df: "pd.DataFrame", cfg: SlalomConfig, labels: List[str]) -> None:
    """Populate df['r'] from the per-population LD columns (or the custom panel)."""
    if cfg.weighted_average_r is not None:
        n_samples = []
        ld = []
        for pop, weight in cfg.weighted_average_r.items():
            if isinstance(weight, str):
                if weight not in df.columns:
                    logger.warning("Column %s not found; skipping in weighted average.", weight)
                    continue
                n_samples.append(df[weight].to_numpy())
            else:
                n_samples.append(np.tile(weight, len(df.index)))
            ld.append(df[f"gnomad_lead_r_{pop}"].to_numpy())
        if not n_samples:
            raise ValueError(
                "No valid --weighted-average-r weights: none of the requested weight columns "
                "were found in the input."
            )
        if len(n_samples) == 1:
            df["r"] = ld[0]
        else:
            weights = np.array(n_samples, dtype=np.float64).T
            ld_arr = np.array(ld, dtype=np.float64).T
            df["r"] = np.nansum(weights * ld_arr, axis=1) / np.nansum(weights * ~np.isnan(ld_arr), axis=1)
    elif cfg.ld_reference == "custom":
        df["r"] = df[labels[0]]
    else:
        df["r"] = df["gnomad_lead_r_nfe"]


def _build_summary(df: "pd.DataFrame", cfg: SlalomConfig) -> "pd.DataFrame":
    """Build the per-locus summary table."""
    import pandas as pd

    df = df.copy()
    df["r2"] = df["r"] ** 2
    n_r2 = int(np.sum(df["r2"] > cfg.r2_threshold))
    n_na = int(np.sum(np.isnan(df["r2"])))
    outlier_idx = (df["r2"] > cfg.r2_threshold) & (df["nlog10p_dentist_s"] > cfg.nlog10p_dentist_s_threshold)
    n_dentist_s_outlier = int(np.sum(outlier_idx))
    max_pip_idx = int(df["prob"].idxmax())
    variant = df["chromosome"].str.cat([df["position"].astype(str), df["allele1"], df["allele2"]], sep=":")
    expr = {
        "lead_pip_variant": [variant.iloc[max_pip_idx]],
        "n_total": [len(df.index)],
        "n_r2": [n_r2],
        "n_na": [n_na],
        "n_dentist_s_outlier": [n_dentist_s_outlier],
        "fraction": [n_dentist_s_outlier / n_r2 if n_r2 > 0 else 0],
        "max_pip": [np.max(df["prob"])],
    }

    if cfg.annotate_consequence:
        nonsyn_idx = (df["r2"] > cfg.r2_threshold) & df["consequence"].isin(resources.NONSYN_CONSEQUENCES)
        expr = {
            **expr,
            "n_nonsyn": [int(np.sum(nonsyn_idx))],
            "n_nonsyn_outlier": [int(np.sum(nonsyn_idx & outlier_idx))],
            "max_pip_nonsyn": [np.max(df["prob"].loc[nonsyn_idx]) if np.any(nonsyn_idx) else np.nan],
            "cs_nonsyn": [bool(np.any(df["cs"].loc[nonsyn_idx]))],
            "cs_99_nonsyn": [bool(np.any(df["cs_99"].loc[nonsyn_idx]))],
            "nonsyn_variants": [",".join(variant.loc[nonsyn_idx].values)],
        }

    if "n_samples" in df.columns:
        if cfg.case_control:
            if "n_cases" in df.columns:
                df["n_eff_samples"] = (
                    df["n_samples"] * (df["n_cases"] / df["n_samples"]) * (1 - df["n_cases"] / df["n_samples"])
                )
            else:
                df["n_eff_samples"] = np.nan
        else:
            df["n_eff_samples"] = df["n_samples"]

        n_eff_r2 = df["n_eff_samples"].loc[df["r2"] > cfg.r2_threshold]
        expr = {
            **expr,
            "min_neff_r2": [np.nanmin(n_eff_r2) if n_r2 > 0 else np.nan],
            "max_neff_r2": [np.nanmax(n_eff_r2) if n_r2 > 0 else np.nan],
        }

    return pd.DataFrame(expr)
