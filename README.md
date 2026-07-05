# SLALOM

SLALOM (<ins>s</ins>uspicious <ins>l</ins>oci <ins>a</ins>na<ins>l</ins>ysis <ins>o</ins>f <ins>m</ins>eta-analysis summary statistics) is a summary statistics-based QC method that identifies suspicious loci for meta-analysis fine-mapping by detecting association statistics outliers based on local LD structure. SLALOM only takes GWAS summary statistics and ancestry-matched external LD reference (e.g., [gnomAD](https://gnomad.broadinstitute.org/downloads#v2-linkage-disequilibrium)) as input and predicts whether each locus shows a suspicious pattern that called into question fine-mapping accuracy. The outlier detection was built upon the simplified version of [the DENTIST method](https://doi.org/10.1038/s41467-021-27438-7).

Analysis and figure generation code for [Kanai, M. et al. (2022)](http://dx.doi.org/10.1101/2022.03.16.22272457) is available [here](https://github.com/mkanai/slalom-paper). Fine-mapping pipeline is available [here](https://github.com/mkanai/finemapping-pipeline).

<p align="center"><img src="https://mkanai.github.io/assets/img/slalom.svg" width="40%"></p>

## What's new: Hail-free

SLALOM is now a pip-installable package with a `slalom` command-line tool and **needs no
Hail or Spark at runtime**. The gnomAD LD `BlockMatrix` is read directly from cloud storage
in pure Python via [`ldcov`](https://github.com/mkanai/ldcov), and every annotation that
used to read a Hail `.ht` (allele alignment, CUPs, VEP consequence, gnomAD frequency) now
reads a Parquet reference table instead. See [`CHANGELOG.md`](CHANGELOG.md) for the full
list of changes from the original single-file `slalom.py`.

## Installation

```bash
pip install slalom-qc
```

The distribution is published on PyPI as `slalom-qc`; the import package and CLI are both
`slalom` (`import slalom`, `slalom --help`).

Python ≥ 3.9. Reading LD or reference tables from Google Cloud Storage uses `gcsfs`
(installed automatically); reading a `BlockMatrix` from AWS S3 (e.g. Pan-UKB) needs the
`s3` extra: `pip install "slalom-qc[s3]"`.

## Usage

```bash
slalom \
    --snp example/example.snp \
    --out example/example.slalom.txt \
    --out-summary example/example.summary.txt \
    --reference-genome GRCh38 \
    --align-alleles \
    --annotate-consequence \
    --annotate-cups \
    --annotate-gnomad-freq \
    --export-r \
    --lead-variant-choice prob \
    --weighted-average-r afr=n_afr amr=n_amr eas=n_eas fin=n_fin nfe=n_nfe \
    --dentist-s \
    --abf \
    --summary \
    --case-control
```

Run `slalom --help` for the full flag list. The tool can also be used as a library:

```python
from slalom import SlalomConfig, run_slalom

run_slalom(SlalomConfig(snp="example/example.snp", out="out.txt",
                        reference_genome="GRCh38", dentist_s=True))
```

### Reading requester-pays / private buckets

gnomAD LD variant indices and the SLALOM reference tables live in requester-pays buckets.
Pass fsspec options as JSON via `--storage-options`, e.g. to bill your own project:

```bash
slalom ... --storage-options '{"requester_pays": true, "project": "YOUR_GCP_PROJECT"}'
```

## Reference data

All reference paths default to hosted copies and are overridable on the command line.

| Data | Default location | Override flag |
| --- | --- | --- |
| gnomAD LD `BlockMatrix` (r) | `gs://gcp-public-data--gnomad/.../{pop}.common.ld.bm` (public) | n/a |
| gnomAD LD variant index (Parquet) | `gs://ldcov-requester-pays/gnomad_v2.{pop}.{build}.variant_index.parquet` | `--ld-variant-index-dir` |
| gnomAD sites annotations (Parquet) | `gs://finucane-requester-pays/slalom/parquet/...` | `--gnomad-sites-parquet` |
| Conversion-unstable positions (Parquet) | `gs://finucane-requester-pays/slalom/parquet/...` | `--cup-parquet` |

### Downloadable bundles (no requester-pays setup)

The requester-pays reference tables are also published as tar.gz bundles. Download them once
on a machine that has `gcloud` and a billing project, extract, and point SLALOM at the local
copies with the override flags; then runs need no requester-pays billing.

| Bundle | Contents | Location |
| --- | --- | --- |
| `slalom.{build}.parquet.tar.gz` | gnomAD sites + CUP Parquet (b37 ≈ 2.5 GB, b38 ≈ 7.3 GB) | `gs://finucane-requester-pays/slalom/bundles/` |
| `gnomad_v2.{build}.variant_index.tar.gz` | gnomAD LD variant indices (ldcov) | `gs://ldcov-requester-pays/bundles/` |

`{build}` is `b37` (GRCh37) or `b38` (GRCh38); ldcov also publishes Pan-UKB variant-index
bundles (`panukb.{build}.variant_index.tar.gz`).

```bash
BUILD=b38   # or b37
gcloud storage cp \
    gs://finucane-requester-pays/slalom/bundles/slalom.${BUILD}.parquet.tar.gz \
    gs://ldcov-requester-pays/bundles/gnomad_v2.${BUILD}.variant_index.tar.gz \
    . --billing-project YOUR_PROJECT
tar -xzf slalom.${BUILD}.parquet.tar.gz                          # -> sites + CUP directories
mkdir -p ld_index && tar -xzf gnomad_v2.${BUILD}.variant_index.tar.gz -C ld_index

# then point SLALOM at the local paths (GRCh38 shown):
slalom ... --reference-genome GRCh38 \
    --gnomad-sites-parquet gnomad.genomes.r3.1.2.sites.most_severe.b38.parquet \
    --cup-parquet FASTA_BED.ALL_GRCh38.cups.parquet \
    --ld-variant-index-dir ld_index
```

This avoids requester-pays reads, but is **not** a fully offline mode: the gnomAD LD
`BlockMatrix` is still read at runtime from the **public** `gs://gcp-public-data--gnomad`
bucket (free, but GCS network access is required; the matrices are far too large to
download).

The Parquet variant indices are provided by ldcov (see its README for pre-computed bundles
and the `gnomAD`/`Pan-UKB` populations available). The gnomAD sites and CUP Parquet tables
are built **once** from the original Hail Tables with the helpers under `scripts/` (Hail is
only needed for this one-time step; install it with `pip install "slalom-qc[convert]"`):

```bash
# gnomAD sites: written via Spark as a *partitioned Parquet directory* (the genome-wide HT
# has ~760M variants). pyarrow.dataset reads the directory transparently.
python scripts/make_gnomad_sites_parquet.py \
    --ht gs://.../gnomad.genomes.r3.1.2.sites.most_severe.ht \
    --reference-genome GRCh38 --coalesce 256 \
    --out gs://YOUR_BUCKET/.../gnomad.genomes.r3.1.2.sites.most_severe.b38.parquet

python scripts/make_cup_parquet.py \
    --cup-ht    gs://.../FASTA_BED.ALL_GRCh38.novel_CUPs.ht \
    --reject-ht gs://.../FASTA_BED.ALL_GRCh38.reject_2.ht \
    --out FASTA_BED.ALL_GRCh38.cups.parquet
```

### Custom LD reference

To use your own LD instead of gnomAD, point SLALOM at a Hail `BlockMatrix` and its Parquet
variant index (e.g. one built with `ldcov`'s `scripts/make_bm_variant_index.py`):

```bash
slalom ... \
    --ld-reference custom \
    --custom-ld-path my_ld.bm \
    --custom-ld-variant-index-path my_ld.variant_index.parquet \
    --custom-ld-label my_panel
```

## Input format

### (Required) per-locus .snp file

Required minimum columns are as follows:

- `chromosome`: chromosome either in GRCh37 (1, 2, 3...) or in GRCh38 (chr1, chr2, chr3, ...). Users can specify a reference genome by `--reference-genome GRCh38`.
- `position`: position
- `allele1`: reference allele in a specified reference genome. If users are unsure about reference/alternative alleles, set `--align-alleles` to make it consistent with gnomAD.
- `allele2`: alternative allele in a specified reference genome. This allele is assumed to be an effect allele regardless of `--align-alleles`.
- `beta`: effect size
- `se`: standard error
- `p`: P-value

Other input column specifications are as follows:

- If `--weighted-average-r` is specified, sample size columns supplied by this argument are also required, such as `n_afr`, `n_eas`, `n_nfe`, ... (`--weighted-average-r` requires `--export-r`).
- If a total sample size `n_samples` (and `n_cases` for `--case-control`) exist in the input, additional output columns `min_neff_r2` and `max_neff_r2` will be added.
- Any other input columns will remain in an output except for those overwritten by SLALOM.
- `--summary` reports DENTIST-S outliers and max-PIP statistics, so it requires `--dentist-s` and `--abf`.
- LD variants are matched to the reference panel in exact ref/alt orientation. `--align-alleles` additionally allows matching variants stored in the panel with ref/alt swapped (sign-flipping their `r`), consistent with its role of reconciling allele orientation against gnomAD.

To make SLALOM-compatible per-locus .snp files from a genome-wide summary statistics, you can also use [make_finemap_inputs.py](https://github.com/FINNGEN/finemapping-pipeline/blob/master/python/make_finemap_inputs.py) from our fine-mapping pipeline.

## Legacy WDL pipeline

The `legacy/wdl/` pipeline predates this packaging and targets the standalone Hail/Spark
script. It is kept for reference only; the recommended entry point is now the `slalom` CLI
above.

## Citation

Kanai, M. et al. [Meta-analysis fine-mapping is often miscalibrated at single-variant resolution](http://dx.doi.org/10.1016/j.xgen.2022.100210). Cell Genomics 2, 100210 (2022)

## Contact

Masahiro Kanai (mkanai@broadinstitute.org)
