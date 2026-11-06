"""Node ages from the TimeTree REST API (api.timetree.org), keyed by NCBI taxon id."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterable, Sequence

from .webapi import ApiError, JsonCache, request_json

BASE = "https://api.timetree.org"


class TimeTree:
    def __init__(self, cache: JsonCache | None = None, fetch: Callable = request_json):
        self.cache = cache or JsonCache(None)
        self.fetch = fetch

    @staticmethod
    def _key(ids: list[int]) -> str:
        return "mrca:" + "+".join(map(str, ids))

    def _fetch_one(self, ids: list[int]) -> dict:
        try:
            result = self.fetch(f"{BASE}/mrca/id/{'+'.join(map(str, ids))}/json")
        except ApiError:
            result = None
        found, age = [], None
        if isinstance(result, dict):
            found = [int(i) for i in result.get("found_ids") or []]
            value = result.get("precomputed_age")
            age = float(value) if value not in (None, "", 0) else None
        return {"found": found, "age": age}

    def _query(self, ids: list[int]) -> dict:
        key = self._key(ids)
        if key not in self.cache:
            self.cache.set(key, self._fetch_one(ids))
        return self.cache.get(key)

    def prefetch(self, id_sets: Iterable[Iterable[int]], workers: int = 4) -> None:
        """Fetch many MRCA queries concurrently into the cache (a few at a time, to be polite)."""
        todo = {}
        for ids in id_sets:
            ids = sorted({int(i) for i in ids})
            if len(ids) >= 2 and self._key(ids) not in self.cache:
                todo[self._key(ids)] = ids
        if not todo:
            return
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for key, result in zip(todo, pool.map(self._fetch_one, todo.values())):
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
