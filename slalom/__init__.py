"""SLALOM: suspicious loci analysis of meta-analysis summary statistics.

SLALOM is a summary-statistics-based QC method that flags suspicious loci for
meta-analysis fine-mapping by detecting association-statistic outliers against a
local LD reference (e.g. gnomAD), using a simplified form of DENTIST (DENTIST-S).

This package reads Hail ``BlockMatrix`` LD stores in pure Python via ``ldcov`` and
needs no Hail/Spark at runtime. See :func:`slalom.pipeline.run_slalom`.
"""

from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
from importlib.metadata import version as _version

try:
    __version__ = _version("slalom")
except _PackageNotFoundError:  # running from a source tree without an install
    __version__ = "0.0.0+unknown"

# Lazy attribute loading keeps `import slalom` cheap (numpy/pandas/pyarrow are only
# pulled in when a symbol that needs them is first accessed).
_LAZY_IMPORTS = {
    "run_slalom": ".pipeline",
    "SlalomConfig": ".pipeline",
    "abf": ".stats.abf",
    "get_cs": ".stats.abf",
    "dentist_s": ".stats.dentist",
    "align_alleles": ".annotate.alleles",
    "lead_variant_r": ".ld",
}


def __getattr__(name):
    if name in _LAZY_IMPORTS:
        import importlib

        module = importlib.import_module(_LAZY_IMPORTS[name], package=__name__)
        attr = getattr(module, name)
        globals()[name] = attr  # cache for subsequent access
        return attr
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
