"""Region-filtered readers for the Parquet reference tables (gnomAD sites, CUPs).

Each table is a Parquet dataset keyed by ``contig``/``position``; a per-locus query pulls
only the relevant slice via pyarrow predicate pushdown, so no Hail/Spark is needed at
runtime. Build the Parquet copies once with the helpers under ``scripts/`` (see README).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional

import pyarrow.compute as pc
import pyarrow.dataset as ds
from ldcov.io.fs_utils import resolve_filesystem

if TYPE_CHECKING:
    import pandas as pd


class _ParquetTable:
    """Thin wrapper over a Parquet dataset with a contig/position region query."""

    def __init__(self, path: str, storage_options: Optional[dict] = None) -> None:
        self.path = path
        fs, inner = resolve_filesystem(path, storage_options)
        self._dataset = ds.dataset(inner, format="parquet", filesystem=fs)

    @property
    def columns(self) -> List[str]:
        return self._dataset.schema.names

    def query_region(self, chrom: str, start: int, end: int, columns: Optional[List[str]] = None) -> "pd.DataFrame":
        flt = (
            (pc.field("contig") == str(chrom))
            & (pc.field("position") >= int(start))
            & (pc.field("position") <= int(end))
        )
        table = self._dataset.to_table(columns=columns, filter=flt)
        return table.to_pandas()


class ReferencePanel:
    """Lazily-opened Parquet reference tables for a single reference genome build.

    Only the tables actually needed by the requested annotations are opened, so a run
    that skips (say) CUP annotation never touches the CUP Parquet.
    """

    def __init__(
        self,
        sites_path: Optional[str] = None,
        cup_path: Optional[str] = None,
        storage_options: Optional[dict] = None,
    ) -> None:
        self._sites_path = sites_path
        self._cup_path = cup_path
        self._storage_options = storage_options
        self._sites: Optional[_ParquetTable] = None
        self._cups: Optional[_ParquetTable] = None

    @property
    def sites(self) -> _ParquetTable:
        if self._sites is None:
            if self._sites_path is None:
                raise ValueError("gnomAD sites Parquet path is not configured")
            self._sites = _ParquetTable(self._sites_path, self._storage_options)
        return self._sites

    @property
    def cups(self) -> _ParquetTable:
        if self._cups is None:
            if self._cup_path is None:
                raise ValueError("CUP Parquet path is not configured")
            self._cups = _ParquetTable(self._cup_path, self._storage_options)
        return self._cups

    def query_sites(self, chrom: str, start: int, end: int, columns: Optional[List[str]] = None) -> "pd.DataFrame":
        """gnomAD sites annotation rows on `chrom` with start <= position <= end."""
        return self.sites.query_region(chrom, start, end, columns=columns)

    def query_cups(self, chrom: str, start: int, end: int) -> "pd.DataFrame":
        """CUP intervals on `chrom` that overlap the half-open window [start, end].

        The CUP table stores half-open intervals [start, end); an interval overlaps the
        query window when its start <= end-of-window and its end > start-of-window.
        """
        flt = (pc.field("contig") == str(chrom)) & (pc.field("start") <= int(end)) & (pc.field("end") > int(start))
        table = self.cups._dataset.to_table(columns=["contig", "start", "end"], filter=flt)
        return table.to_pandas()
