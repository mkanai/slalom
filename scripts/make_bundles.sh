#!/usr/bin/env bash
# Package the SLALOM reference Parquet tables into per-build tar.gz bundles for offline use.
#
# Each bundle holds the gnomAD sites + CUP Parquet directories for one reference genome, so a
# user can download it on a machine with GCS billing configured, extract it, and copy the
# directories to a cluster without requester-pays access. Point --gnomad-sites-parquet /
# --cup-parquet at the extracted directories there.
#
# Usage:
#   scripts/make_bundles.sh <src_parquet_prefix> <dst_bundle_prefix> <billing_project> [workdir]
# Example:
#   scripts/make_bundles.sh \
#       gs://finucane-requester-pays/slalom/parquet \
#       gs://finucane-requester-pays/slalom/bundles \
#       my-gcp-project
set -euo pipefail

SRC="${1:?src parquet prefix}"; SRC="${SRC%/}"
DST="${2:?dst bundle prefix}"; DST="${DST%/}"
PROJECT="${3:?billing project}"
WORK="${4:-$(mktemp -d)}"
GZIP="$(command -v pigz || echo gzip)"

# build tag -> (sites dir, cups dir)
declare -A SITES=(
  [b37]="gnomad.genomes.r2.1.1.sites.most_severe.b37.parquet"
  [b38]="gnomad.genomes.r3.1.2.sites.most_severe.b38.parquet"
)
declare -A CUPS=(
  [b37]="FASTA_BED.ALL_GRCh37.cups.parquet"
  [b38]="FASTA_BED.ALL_GRCh38.cups.parquet"
)

for build in b37 b38; do
  echo "[$(date -u +%H:%M:%S)] building $build bundle in $WORK"
  stage="$WORK/$build"; mkdir -p "$stage"
  gcloud storage cp -r "$SRC/${SITES[$build]}" "$stage/" --billing-project="$PROJECT"
  gcloud storage cp -r "$SRC/${CUPS[$build]}"  "$stage/" --billing-project="$PROJECT"
  tarball="$WORK/slalom.$build.parquet.tar.gz"
  tar -C "$stage" -cf - "${SITES[$build]}" "${CUPS[$build]}" | "$GZIP" > "$tarball"
  echo "[$(date -u +%H:%M:%S)] $(du -h "$tarball" | cut -f1) -> $DST/slalom.$build.parquet.tar.gz"
  gcloud storage cp "$tarball" "$DST/slalom.$build.parquet.tar.gz" --billing-project="$PROJECT"
  rm -rf "$stage" "$tarball"
done
echo "[$(date -u +%H:%M:%S)] done"
