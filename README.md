# SLALOM

SLALOM (<ins>s</ins>uspicious <ins>l</ins>oci <ins>a</ins>na<ins>l</ins>ysis <ins>o</ins>f <ins>m</ins>eta-analysis summary statistics) is a summary statistics-based QC method that identifies suspicious loci for meta-analysis fine-mapping by detecting association statistics outliers based on local LD structure. SLALOM only takes GWAS summary statistics and ancestry-matched external LD reference (e.g., [gnomAD](https://gnomad.broadinstitute.org/downloads#v2-linkage-disequilibrium)) as input and predicts whether each locus shows a suspicious pattern that called into question fine-mapping accuracy. The outlier detection was built upon the simplified version of [the DENTIST method](https://doi.org/10.1038/s41467-021-27438-7).

Analysis and figure generation code for [Kanai, M. et al. (2022)](http://dx.doi.org/10.1101/2022.03.16.22272457) is available [here](https://github.com/mkanai/slalom-paper). Fine-mapping pipeline is available [here](https://github.com/mkanai/finemapping-pipeline).

<p align="center"><img src="https://mkanai.github.io/assets/img/slalom.svg" width="40%"></p>

## Requirements

- Python 3.7 or later
- [Hail v0.2](https://hail.is/)
- numpy
- scipy
- pandas

To run our WDL pipeline on Google Cloud, you additionally need:

- [Cromwell](https://cromwell.readthedocs.io/en/stable/)
- Active Google Cloud project
  - Note: A part of reference files are located in a public requester-pays bucket (`gs://finucane-requester-pays`)

To run SLALOM locally, you need:

- [Cloud Storage Connector](https://hail.is/docs/0.2/cloud/google_cloud.html#reading-from-google-cloud-storage)

The following command would be the easiest way of installation.

```bash
curl -sSL https://broad.io/install-gcs-connector | python3 - --gcs-requester-pays-project YOUR_PROJECT_ID
```

## Usage

### (Recommended) WDL pipeline

Please modify `wdl/slalom_example.json` and submit with `wdl/slalom.wdl` and `wdl/slalom_sub.zip`.

### Running per-locus

Example files are available at `./example` which was created from the GBMI meta-analysis summary statistics for COPD available [here](https://www.globalbiobankmeta.org/resources).

```bash
PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=1g pyspark-shell" \
python3 slalom.py \
        --snp example/example.snp \
        --out example/example.slalom.txt \
        --out-summary example/example.summary.txt \
        --annotate-consequence \
        --annotate-cups \
        --annotate-gnomad-freq \
        --export-r \
        --lead-variant-choice "prob" \
        --weighted-average-r afr=n_afr amr=n_amr eas=n_eas fin=n_fin nfe=n_nfe \
        --dentist-s \
        --abf \
        --summary \
        --case-control \
        --reference-genome GRCh38
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

- If `--weighted-average-r` is specified, sample size columns supplied by this argument are also required, such as `n_afr`, `n_eas`, `n_nfe`, ...
- If a total sample size `n_samples` (and `n_cases` for `--case-control`) exist in the input, additional output columns `min_neff_r2` and `max_neff_r2` will be added.
- Any other input columns will remain in an output except for those overwritten by SLALOM.

To make SLALOM-compatible per-locus .snp files from a genome-wide summary statistics, you can also use [make_finemap_inputs.py](https://github.com/FINNGEN/finemapping-pipeline/blob/master/python/make_finemap_inputs.py) from our fine-mapping pipeline.

### (Optional) WDL input (.json)

To use our WDL pipeline, please modify `wdl/slalom_example.json`. Specifications for the following options are as follows:

- `slalom.sumstats_pattern`: Path pattern for a summary statistics where `{PHENO}` will be repalced by a phenotype name. E.g., `gs://YOUR_BUCKET/{PHENO}.sumstats.txt.gz`.
- `slalom.phenolistfile`: Path to a plain text file without header. The first column corresponds to a phenotype name (`{PHENO}` above). The second column corresponds to an argument for ``--weighted-average-r`.

## Citation

Kanai, M. et al. [Meta-analysis fine-mapping is often miscalibrated at single-variant resolution](http://dx.doi.org/10.1016/j.xgen.2022.100210). Cell Genomics 2, 100210 (2022)

## Contact

Masahiro Kanai (mkanai@broadinstitute.org)
