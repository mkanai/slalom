"""Statistical primitives for SLALOM (ABF fine-mapping and the DENTIST-S outlier test)."""

from .abf import abf, get_cs
from .dentist import dentist_s

__all__ = ["abf", "get_cs", "dentist_s"]
