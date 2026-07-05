"""Filesystem-agnostic IO for SLALOM (local paths and gs://, backed by fsspec)."""

from .reference import ReferencePanel
from .snp import read_snp, write_table

__all__ = ["read_snp", "write_table", "ReferencePanel"]
