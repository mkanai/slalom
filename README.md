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

To run SLALOM locally, you need:
- [Cloud Storage Connector](https://hail.is/docs/0.2/cloud/google_cloud.html#reading-from-google-cloud-storage)

The following command would be the easiest way of installation.
```bash
curl -sSL https://broad.io/install-gcs-connector | python3
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

## Citation

Kanai, M. et al. [Meta-analysis fine-mapping is often miscalibrated at single-variant resolution](http://dx.doi.org/10.1101/2022.03.16.22272457). medRxiv (2022) doi:10.1101/2022.03.16.22272457

## Contact

Masahiro Kanai (mkanai@broadinstitute.org)
