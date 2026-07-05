"""Hail-free variant annotations backed by Parquet reference tables."""

from .alleles import align_alleles
from .annotations import annotate_consequence_and_freq, annotate_cups

__all__ = ["align_alleles", "annotate_cups", "annotate_consequence_and_freq"]
