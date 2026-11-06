"""Node ages from the TimeTree REST API (api.timetree.org), keyed by NCBI taxon id."""

from __future__ import annotations

from typing import Callable, Sequence

from .webapi import ApiError, JsonCache, request_json

BASE = "https://api.timetree.org"


class TimeTree:
    def __init__(self, cache: JsonCache | None = None, fetch: Callable = request_json):
        self.cache = cache or JsonCache(None)
        self.fetch = fetch

    def _query(self, ids: list[int]) -> dict:
        key = "mrca:" + "+".join(map(str, ids))
        if key not in self.cache:
            try:
                result = self.fetch(f"{BASE}/mrca/id/{'+'.join(map(str, ids))}/json")
            except ApiError:
                result = None
            found, age = [], None
            if isinstance(result, dict):
                found = [int(i) for i in result.get("found_ids") or []]
                value = result.get("precomputed_age")
                age = float(value) if value not in (None, "", 0) else None
            self.cache.set(key, {"found": found, "age": age})
        return self.cache.get(key)

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
