# Changelog

All notable changes to this project are documented in this file. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [2.0.0] - 2026-07-04

First packaged, Hail-free release on PyPI (major bump from the standalone `slalom.py` script,
which was distributed as 1.0.0). SLALOM is now a pip-installable package with a `slalom`
console entry point.

### Added
- Modern package layout (`slalom.{stats,io,annotate,ld,cli,pipeline}`) with a `slalom` CLI.
- Programmatic API: `slalom.run_slalom(SlalomConfig(...))`.
- One-time Hail conversion helpers under `scripts/` to build the Parquet reference tables,
  and `scripts/make_bundles.sh` to package them into downloadable tar.gz bundles.
- `--ld-variant-index-dir` to read gnomAD LD variant indices from a local directory (e.g. an
  extracted ldcov bundle), enabling runs with no requester-pays access.
- Inline type hints across the package and a PEP 561 `py.typed` marker so downstream
  type-checkers pick up the annotations, plus a lenient `mypy` config and a CI `typecheck` job.

### Changed
- Run invariants (custom-LD paths; `weighted_average_r` requires `export_r`; `summary`
  requires `dentist_s` and `abf`) are now validated in `SlalomConfig` construction, so the
  programmatic `run_slalom(SlalomConfig(...))` API rejects invalid configs the same way the
  CLI does. Previously these were only checked in the CLI layer.
- **LD is read from the gnomAD Hail `BlockMatrix` in pure Python via
  [`ldcov`](https://github.com/mkanai/ldcov)** — no Hail or Spark at runtime. Only the
  blocks along the lead variant's row/column are read.
- Variant → matrix-index mapping uses ldcov's hosted Parquet variant indices.
- Allele alignment, CUP, consequence, and gnomAD-frequency annotations were reimplemented
  against Parquet reference tables (region-filtered pyarrow queries + pandas), replacing
  the Hail-table joins. Output column layout is preserved.
- Input `.snp` reading and output writing use pandas + fsspec (local and `gs://`).

### Fixed
- Default `export_r=False` runs no longer crash: LD columns are kept as signed `r`
  internally (so DENTIST-S sees signed `r`, not `r^2`) and converted to `r^2` output columns
  only at the end. Previously the gnomAD default raised `KeyError` and the custom path fed
  `r^2` into DENTIST-S (a latent bug inherited from the script, which was only run with
  `--export-r`).
- `get_cs` returns an empty credible set instead of raising when PIPs are all non-finite.
- `--weighted-average-r` raises a clear error when none of the requested weight columns exist.

### Removed
- The Hail dependency and the standalone `slalom.py` script.
- The unused `--delimiter` flag.

### Notes
- `--summary` now explicitly requires `--abf` and `--dentist-s` (it always consumed their
  outputs); `--weighted-average-r` requires `--export-r`.
- LD variant matching is delegated to ldcov's `VariantIndex.match_variants` (>= ldcov
  0.5.0), using exact ref/alt orientation. `--align-alleles` additionally enables ldcov's
  ref/alt-swap matching (sign-flipping r), consistent with its role of reconciling allele
  orientation against gnomAD.
- Validated against the bundled `example/`: every computed column (LD, DENTIST-S, PIP,
  credible sets, consequence, CUPs, summary) is bit-identical to the original output. The
  gnomAD `af_*` annotation columns now carry full float64 precision instead of the original
  Hail TSV export's ~5-significant-figure rounding (max abs diff < 5e-6; no downstream
  effect).

## [1.0.0] - 2022

Standalone Hail-based `slalom.py` script accompanying Kanai et al. (2022).
