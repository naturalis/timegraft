"""Node ages from the TimeTree REST API (api.timetree.org), keyed by NCBI taxon id.

Cache format
------------
Each MRCA query is cached under ``mrca:<id>+<id>+...``. Version 1 entries
hold only what the dating rules use, ``{"found": [...], "age": ...}``.
Version 2 entries (``"v": 2``) keep the rest of TimeTree's answer as
evidence for later weighting and diagnostics (#47): the individual time
estimates, the confidence interval, the studies behind them, and the NCBI
id of the MRCA TimeTree dated. ``found`` and ``age`` mean the same in both
versions, and the dating rules read only those two.

A version 1 entry can be upgraded in place (``upgrade_legacy``): TimeTree is
asked again for the evidence, but the cached ``found`` and ``age`` are kept,
so that ages do not change with the cache format. If the fresh answer
differs, it is recorded in ``fresh_found`` and ``fresh_age``.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterable, Sequence

from .webapi import ApiError, JsonCache, request_json

BASE = "https://api.timetree.org"
CACHE_VERSION = 2

# Columns of :meth:`TimeTree.evidence`, in the order they are written to the outputs.
EVIDENCE_COLUMNS = (
    "evidence", "found_ids", "mrca_ncbi", "n_estimates", "n_studies", "estimate_min", "estimate_max",
    "ci_low", "ci_high", "api_version",
)

log = logging.getLogger(__name__)


def _float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _int(value) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _estimates(value) -> list[float]:
    """TimeTree's ``time_estimates``: a Postgres array literal such as ``"{86,91,94.5}"``."""
    if isinstance(value, (list, tuple)):
        items = value
    elif isinstance(value, str):
        items = value.strip().strip("{}").split(",")
    else:
        return []
    return [f for f in (_float(v) for v in items) if f is not None]


def parse_response(result) -> dict:
    """A cache entry (version 2) from one answer of ``/mrca/id/.../json``.

    ``found`` and ``age`` follow the version 1 rules exactly: ``age`` is
    ``precomputed_age`` unless that is missing or zero.
    """
    entry: dict = {"v": CACHE_VERSION, "found": [], "age": None}
    if not isinstance(result, dict):
        entry["error"] = "no answer"
        return entry
    entry["found"] = [int(i) for i in result.get("found_ids") or []]
    value = result.get("precomputed_age")
    entry["age"] = float(value) if value not in (None, "", 0) else None
    entry["missing"] = [int(i) for i in result.get("missing_ids") or []]
    entry["mrca_ncbi"] = _int(result.get("taxon_id"))
    entry["mrca_ttid"] = _int(result.get("mrca_ttid"))
    entry["estimates"] = _estimates(result.get("time_estimates"))
    low, high = _float(result.get("precomputed_ci_low")), _float(result.get("precomputed_ci_high"))
    # TimeTree reports 0, 0 when it has no interval (a single estimate)
    entry["ci_low"], entry["ci_high"] = (low, high) if (low or high) else (None, None)
    entry["adjusted_age"] = _float(result.get("adjusted_age")) or None
    entry["studies"] = [
        {
            "author": s.get("first_author"),
            "year": _int(s.get("pub_year")),
            "trees": _int(s.get("tree_count")),
            "age": _float(s.get("age_estimate")),
        }
        for s in result.get("study_data") or []
        if isinstance(s, dict)
    ]
    if result.get("error"):
        entry["error"] = str(result["error"])
    if result.get("version") is not None:
        entry["api_version"] = str(result["version"])
    return entry


def is_legacy(entry) -> bool:
    return isinstance(entry, dict) and entry.get("v", 1) < CACHE_VERSION


class TimeTree:
    def __init__(self, cache: JsonCache | None = None, fetch: Callable = request_json, upgrade_legacy: bool = False):
        self.cache = cache or JsonCache(None)
        self.fetch = fetch
        self.upgrade_legacy = upgrade_legacy

    @staticmethod
    def _key(ids: list[int]) -> str:
        return "mrca:" + "+".join(map(str, ids))

    def _fetch_one(self, ids: list[int]) -> dict:
        try:
            result = self.fetch(f"{BASE}/mrca/id/{'+'.join(map(str, ids))}/json")
        except ApiError as err:
            entry = parse_response(None)
            entry["error"] = f"HTTP {err.status}"
            return entry
        return parse_response(result)

    def _upgrade(self, old: dict, ids: list[int]) -> dict:
        """A version 2 entry for a version 1 one, keeping its ``found`` and ``age``."""
        try:
            new = self._fetch_one(ids)
        except Exception as err:  # network trouble: keep the old entry, try again next run
            log.warning("could not upgrade the TimeTree cache entry for %s: %s", self._key(ids), err)
            return old
        fresh_found, fresh_age = new["found"], new["age"]
        new["found"], new["age"] = list(old.get("found") or []), old.get("age")
        new["upgraded_from"] = old.get("v", 1)
        if sorted(fresh_found) != sorted(new["found"]):
            new["fresh_found"] = fresh_found
        if fresh_age != new["age"]:
            new["fresh_age"] = fresh_age
        return new

    def _needs_fetch(self, key: str) -> bool:
        return key not in self.cache or (self.upgrade_legacy and is_legacy(self.cache.get(key)))

    def _query(self, ids: list[int]) -> dict:
        ids = sorted({int(i) for i in ids})
        key = self._key(ids)
        if key not in self.cache:
            self.cache.set(key, self._fetch_one(ids))
        elif self.upgrade_legacy and is_legacy(self.cache.get(key)):
            self.cache.set(key, self._upgrade(self.cache.get(key), ids))
        return self.cache.get(key)

    def prefetch(self, id_sets: Iterable[Iterable[int]], workers: int = 4) -> None:
        """Fetch many MRCA queries concurrently into the cache (a few at a time, to be polite)."""
        todo = {}
        for ids in id_sets:
            ids = sorted({int(i) for i in ids})
            if len(ids) >= 2 and self._needs_fetch(self._key(ids)):
                todo[self._key(ids)] = ids
        if not todo:
            return

        def one(ids: list[int]) -> dict:
            old = self.cache.get(self._key(ids))
            return self._upgrade(old, ids) if old is not None else self._fetch_one(ids)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            for key, result in zip(todo, pool.map(one, todo.values())):
                self.cache.set(key, result)

    def split_ids(self, sides: Sequence[Sequence[int]]) -> list[int]:
        return sorted({int(i) for side in sides for i in side})

    def split_age(self, sides: Sequence[Sequence[int]]) -> float | None:
        """Age in Ma of the node that joins the given sides (e.g. two child clades).

        Only accepted when TimeTree found taxa from at least two sides;
        otherwise the MRCA it reports would be a shallower node.
        """
        sides = [sorted({int(i) for i in side}) for side in sides]
        ids = sorted({i for side in sides for i in side})
        if len(ids) < 2:
            return None
        result = self._query(ids)
        found = set(result["found"])
        if sum(1 for side in sides if found & set(side)) < 2:
            return None
        return result["age"]

    def evidence(self, ids: Iterable[int]) -> dict[str, str]:
        """What TimeTree said about a query already made, as strings for a table (see EVIDENCE_COLUMNS).

        Reads the cache only. ``evidence`` is ``full`` for a version 2 entry,
        ``legacy_cache`` for a version 1 entry (only ``found_ids`` is known),
        and ``none`` when the query is not in the cache.
        """
        ids = sorted({int(i) for i in ids})
        entry = self.cache.get(self._key(ids)) if len(ids) >= 2 else None
        out = {c: "" for c in EVIDENCE_COLUMNS}
        if not isinstance(entry, dict):
            out["evidence"] = "none"
            return out
        out["found_ids"] = "+".join(map(str, sorted(entry.get("found") or [])))
        if is_legacy(entry):
            out["evidence"] = "legacy_cache"
            return out
        out["evidence"] = "full"
        estimates = entry.get("estimates") or []

        def fmt(value) -> str:
            return "" if value is None else f"{value:.4f}"

        out.update(
            mrca_ncbi="" if entry.get("mrca_ncbi") is None else str(entry["mrca_ncbi"]),
            n_estimates=str(len(estimates)),
            n_studies=str(len(entry.get("studies") or [])),
            estimate_min=fmt(min(estimates)) if estimates else "",
            estimate_max=fmt(max(estimates)) if estimates else "",
            ci_low=fmt(entry.get("ci_low")),
            ci_high=fmt(entry.get("ci_high")),
            api_version=entry.get("api_version", ""),
        )
        return out
