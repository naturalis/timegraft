"""Catalogue of Life (COL XR) as the hub for identifier joins (#48).

Each genus (and optionally each species) is matched once against a fixed
COL release on ChecklistBank, with the lineage of its reference records as
classification context. From the COL usage, the NCBI and OTT ids of the
same name are found through ChecklistBank's names index (the ``related``
call), and every candidate is checked against COL's classification so that
homonyms are not joined by name alone.

COL is used for identifiers only, never for placement: the topology still
comes from the OTT and synthesis release files. Ambiguity is kept as
ambiguity: a name with several plausible COL usages, or several NCBI or OTT
candidates that the classification cannot separate, is left unresolved and
reported as such.

All answers are cached (in the same JSON cache as TimeTree), keyed by the
release, so a run with a fixed release and cache is reproducible offline.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterable, Mapping
from urllib.parse import urlencode

from .webapi import ApiError, JsonCache, request_json

BASE = "https://api.checklistbank.org"
NCBI_DATASET = 2169  # NCBI Taxonomy in ChecklistBank
OTT_DATASET = 201890  # OTT in ChecklistBank (3.7.2; OTT ids for topology still come from the release files)

# Ranks whose names are compared between COL and a candidate's classification;
# kingdom and above are too broad to tell homonyms apart.
COMPARE_RANKS = ("phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order",
                 "suborder", "infraorder", "superfamily", "family", "subfamily", "tribe")
LINEAGE_RANKS = ("family", "order", "class", "phylum")
GENUS_MATCH_TYPES = ("exact", "variant", "canonical")
# context taxon of the run to COL kingdom
KINGDOMS = {"Metazoa": "Animalia", "Animalia": "Animalia", "Fungi": "Fungi", "Plantae": "Plantae",
            "Viridiplantae": "Plantae", "Chromista": "Chromista", "Protozoa": "Protozoa"}

# Columns of genera.tsv filled from COL (see docs/technical/provenance.md)
GENUS_COLUMNS = (
    "col_id", "col_release", "col_match_type", "col_name", "col_rank", "col_status", "col_accepted_id",
    "col_context", "col_source", "col_sector", "col_family", "col_order", "col_class", "col_phylum",
    "col_alternatives", "col_confidence",
    "col_ncbi_id", "col_ncbi_candidates", "col_ncbi_path", "col_ott_id", "col_ott_candidates", "col_ott_path",
    "col_ncbi_vs_ott", "col_ott_vs_ott",
)
SPECIES_COLUMNS = ("genus", "species", "col_id", "col_match_type", "col_status", "col_accepted_id", "col_source",
                   "col_confidence", "col_ncbi_id", "col_ncbi_candidates", "col_ncbi_path", "ott_ncbi_id",
                   "col_ncbi_vs_ott")

log = logging.getLogger(__name__)


def _compact_usage(u) -> dict:
    if not isinstance(u, dict):
        return {}
    out = {k: u[k] for k in ("id", "name", "rank", "status", "sectorKey") if u.get(k) is not None}
    if isinstance(u.get("accepted"), dict):
        out["accepted"] = {k: u["accepted"][k] for k in ("id", "name") if u["accepted"].get(k) is not None}
    if u.get("classification"):
        out["classification"] = compact_classification(u["classification"])
    return out


def compact_classification(value) -> list[dict]:
    """Only the ranks timegraft compares (and the kingdom), as rank and name."""
    keep = set(COMPARE_RANKS) | {"kingdom"}
    return [{"rank": c.get("rank"), "name": c.get("name")} for c in value or []
            if isinstance(c, dict) and c.get("rank") in keep]


def compact_match(value) -> dict:
    """The parts of a match answer that timegraft reads; keeps the cache small."""
    if not isinstance(value, dict):
        return {}
    out = {"type": value.get("type")}
    if value.get("usage"):
        out["usage"] = _compact_usage(value["usage"])
    if value.get("alternatives"):
        out["alternatives"] = [_compact_usage(a) for a in value["alternatives"]]
    return out


def compact_related(value) -> list[dict]:
    return [{k: u[k] for k in ("id", "rank", "status") if u.get(k) is not None}
            for u in value or [] if isinstance(u, dict)] if isinstance(value, list) else []


class ChecklistBank:
    """The few ChecklistBank calls timegraft needs, against one dataset release, with caching."""

    def __init__(self, release: str | int, cache: JsonCache | None = None, fetch: Callable = request_json):
        self.cache = cache or JsonCache(None)
        self.fetch = fetch
        self.release = str(release)
        self._info: dict | None = None

    def _get(self, key: str, path: str, params: Mapping | None = None, default=None, compact=None):
        """GET ``path`` once; the answer (passed through ``compact``) is cached under ``key``.
        Client errors are cached as ``default``."""
        if key in self.cache:
            return self.cache.get(key)
        url = f"{BASE}{path}" + (("?" + urlencode(params, doseq=True)) if params else "")
        try:
            value = self.fetch(url)
        except ApiError as err:
            if err.status >= 500:
                raise
            value = default
        if compact is not None and value is not None:
            value = compact(value)
        self.cache.set(key, value)
        return value

    @property
    def info(self) -> dict:
        """Key, alias and version of the release. An alias such as ``3LXR`` is resolved to a fixed key."""
        if self._info is None:
            data = self._get(f"col:dataset:{self.release}", f"/dataset/{self.release}", default={},
                             compact=lambda v: {k: v.get(k) for k in ("key", "alias", "title", "version", "issued",
                                                                      "origin")} if isinstance(v, dict) else {}) or {}
            self._info = {k: data.get(k) for k in ("key", "alias", "title", "version", "issued", "origin")}
            if not self._info.get("key"):
                raise ValueError(f"ChecklistBank has no dataset {self.release!r}")
            if str(self._info["key"]) != self.release:
                log.warning("COL release %s resolves to dataset %s (%s); set col.release to %s to keep runs "
                            "reproducible", self.release, self._info["key"], self._info["alias"], self._info["key"])
        return self._info

    @property
    def key(self) -> int:
        return int(self.info["key"])

    def match(self, name: str, rank: str, context: Mapping[str, str]) -> dict:
        params = {"q": name, "rank": rank, "verbose": "true", **{k: v for k, v in sorted(context.items()) if v}}
        cache_key = f"col:{self.key}:match:" + urlencode(params)
        return self._get(cache_key, f"/dataset/{self.key}/match/nameusage", params, default={},
                         compact=compact_match) or {}

    def related(self, usage_id: str, dataset: int) -> list[dict]:
        cache_key = f"col:{self.key}:related:{usage_id}:{dataset}"
        value = self._get(cache_key, f"/dataset/{self.key}/nameusage/{usage_id}/related", {"datasetKey": dataset},
                          default=[], compact=compact_related)
        return value if isinstance(value, list) else []

    def classification(self, dataset: int, usage_id: str) -> list[dict]:
        cache_key = f"col:classification:{dataset}:{usage_id}"
        value = self._get(cache_key, f"/dataset/{dataset}/taxon/{usage_id}/classification", default=[],
                          compact=compact_classification)
        return value if isinstance(value, list) else []

    def sector_mode(self, sector: int | None) -> str:
        """``attach`` or ``union`` for the base release, ``merge`` for sources added in XR."""
        if sector is None:
            return ""
        value = self._get(f"col:{self.key}:sector:{sector}", f"/dataset/{self.key}/sector/{sector}", default={},
                          compact=lambda v: {"mode": v.get("mode")} if isinstance(v, dict) else {})
        return (value or {}).get("mode") or ""


def _names(classification: Iterable[Mapping]) -> set[str]:
    return {c.get("name") for c in classification if c.get("rank") in COMPARE_RANKS and c.get("name")}


def _crosswalk(clb: ChecklistBank, usage_id: str, col_names: set[str], dataset: int, rank: str) -> dict:
    """The id in ``dataset`` of the COL usage, through the names index, checked against COL's classification.

    ``path`` is ``names_index`` (one candidate, classification agrees),
    ``names_index+classification`` (several candidates, one agrees best),
    ``disagrees`` (candidates, none sharing a taxon with COL), ``ambiguous``
    (several agree equally) or ``none`` (the name is not in ``dataset``).
    """
    candidates = [u for u in clb.related(usage_id, dataset) if u.get("rank") == rank]
    accepted = [u for u in candidates if u.get("status") == "accepted"] or candidates
    out = {"id": "", "candidates": ";".join(sorted(str(u["id"]) for u in accepted)), "path": "none"}
    if not accepted:
        return out
    scored = []
    for u in accepted:
        shared = len(col_names & _names(clb.classification(dataset, str(u["id"]))))
        scored.append((shared, str(u["id"])))
    best = max(s for s, _ in scored)
    winners = [i for s, i in scored if s == best]
    if best == 0:
        out["path"] = "disagrees"
    elif len(winners) > 1:
        out["path"] = "ambiguous"
    else:
        out["id"] = winners[0]
        out["path"] = "names_index" if len(accepted) == 1 else "names_index+classification"
    return out


ACCEPTED = ("accepted", "provisionally accepted")


def _pick(result: Mapping, rank: str, context_names: set[str]) -> tuple[dict | None, str, list[str]]:
    """The usage a match result stands for, how it was chosen, and the other candidates.

    ChecklistBank calls a match ``ambiguous`` whenever other usages of the
    name exist, even if they are only synonyms. A usage it does not call
    ambiguous is taken (``single``, or ``matcher`` if other accepted usages
    of the name exist). For an ambiguous one, among the usages of the right
    rank that are accepted (or provisionally so), a single one is taken
    (``single``); of several, the one whose classification shares most taxa
    with ``context_names`` (``context``), if there is exactly one; otherwise
    none. A usage of another rank (``higherrank``) is no match.
    """
    usage = result.get("usage") or {}
    if result.get("type") not in GENUS_MATCH_TYPES + ("ambiguous",):
        return None, "", []
    pool = [u for u in [usage] + list(result.get("alternatives") or []) if u and u.get("rank") == rank]
    seen, unique = set(), []
    for u in pool:
        if u.get("id") not in seen:
            seen.add(u.get("id"))
            unique.append(u)
    accepted = [u for u in unique if u.get("status") in ACCEPTED]
    others = lambda chosen: sorted(str(u.get("id")) for u in accepted if u is not chosen)  # noqa: E731
    if result.get("type") in GENUS_MATCH_TYPES and usage.get("rank") == rank:
        # ChecklistBank's matcher chose; other accepted usages of the name are reported
        return usage, "single" if len(accepted) <= 1 or usage not in accepted else "matcher", others(usage)
    if len(accepted) == 1:
        return accepted[0], "single", []
    if not accepted:
        return None, "", []
    scored = sorted(((len(context_names & _names(u.get("classification") or [])), i) for i, u in enumerate(accepted)),
                    reverse=True)
    if scored[0][0] > 0 and (len(scored) == 1 or scored[0][0] > scored[1][0]):
        chosen = accepted[scored[0][1]]
        return chosen, "context", others(chosen)
    return None, "", others(None)


def _resolve_usage(clb: ChecklistBank, name: str, rank: str, context: dict[str, str], kingdom: str | None) -> dict:
    """Match a name in COL, with the full context first and the kingdom alone if that fails."""
    attempts = [("records", {**context, "kingdom": kingdom or ""})] if any(context.values()) else []
    attempts.append(("kingdom", {"kingdom": kingdom or ""}) if kingdom else ("none", {}))
    context_names = {v for v in context.values() if v}
    found = {"result": {}, "usage": None, "how": "", "alternatives": [], "context": "none"}
    for label, ctx in attempts:
        result = clb.match(name, rank, ctx)
        usage, how, alternatives = _pick(result, rank, context_names)
        found = {"result": result, "usage": usage, "how": how, "alternatives": alternatives, "context": label}
        if usage is not None:
            break
    return found


def _confidence(status: str, context: str, source: str, crosswalk_path: str, how: str) -> str:
    """``high``, ``medium`` or ``low``. See docs/technical/provenance.md."""
    if source == "xr":
        return "low"
    if (status not in ACCEPTED or context != "records" or how in ("context", "matcher")
            or crosswalk_path == "names_index+classification"):
        return "medium"
    return "high"


def resolve_name(clb: ChecklistBank, name: str, rank: str, context: dict[str, str], kingdom: str | None,
                 ott: bool = True) -> dict:
    """One genus or species resolved in COL and crosswalked to NCBI (and OTT), as genera.tsv columns."""
    row = {c: "" for c in GENUS_COLUMNS}
    row["col_release"] = str(clb.key)
    found = _resolve_usage(clb, name, rank, context, kingdom)
    result = found["result"]
    match_type = result.get("type") or "none"
    usage = found["usage"]
    row["col_match_type"] = match_type
    row["col_context"] = found["context"]
    row["col_alternatives"] = ";".join(found["alternatives"])
    if usage is None:
        # nothing, only a higher rank, or several usages that the context cannot separate: unresolved
        returned = result.get("usage") or {}
        if match_type in GENUS_MATCH_TYPES and returned.get("rank") != rank:
            row["col_match_type"] = "higherrank"
        return row
    accepted = usage.get("accepted") or {}
    target = accepted if usage.get("status") != "accepted" and accepted.get("id") else usage
    classification = usage.get("classification") or []
    by_rank = {c.get("rank"): c.get("name") for c in classification}
    mode = clb.sector_mode(usage.get("sectorKey"))
    row.update(
        col_id=str(target.get("id") or usage.get("id")),
        col_name=usage.get("name") or "",
        col_rank=usage.get("rank") or "",
        col_status=usage.get("status") or "",
        col_accepted_id=str(accepted.get("id") or "") if usage.get("status") != "accepted" else "",
        col_source="xr" if mode == "merge" else "base",
        col_sector=str(usage.get("sectorKey") or ""),
        **{f"col_{r}": by_rank.get(r) or "" for r in LINEAGE_RANKS},
    )
    col_names = _names(classification)
    ncbi = _crosswalk(clb, row["col_id"], col_names, NCBI_DATASET, rank)
    row.update(col_ncbi_id=ncbi["id"], col_ncbi_candidates=ncbi["candidates"], col_ncbi_path=ncbi["path"])
    if ott:
        ottx = _crosswalk(clb, row["col_id"], col_names, OTT_DATASET, rank)
        row.update(col_ott_id=ottx["id"], col_ott_candidates=ottx["candidates"], col_ott_path=ottx["path"])
    row["col_confidence"] = _confidence(row["col_status"], row["col_context"], row["col_source"], ncbi["path"],
                                        found["how"])
    return row


def compare_ids(col_id: str, ott_ids: str, first: str) -> str:
    """How a COL-derived id relates to OTT's: ``same`` as the id OTT uses, ``merged`` (another id OTT holds
    for the taxon), ``different``, or empty if either side has none."""
    if not col_id or not (ott_ids or first):
        return ""
    if col_id == str(first or ""):
        return "same"
    if col_id in set(filter(None, (ott_ids or "").split(";"))):
        return "merged"
    return "different"


def resolve_genera(clb: ChecklistBank, genera: Iterable[str], lineage: Mapping[str, Mapping[str, str]],
                   kingdom: str | None, workers: int = 4) -> dict[str, dict]:
    """genera.tsv columns from COL for every genus, with ``lineage`` (records_<rank>) as context."""
    genera = sorted(set(genera))
    clb.info  # resolve the release once, before the threads start

    def one(genus: str) -> dict:
        ctx = {r: (lineage.get(genus) or {}).get(f"records_{r}", "") for r in LINEAGE_RANKS}
        return resolve_name(clb, genus, "genus", ctx, kingdom)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        return dict(zip(genera, pool.map(one, genera)))


def resolve_species(clb: ChecklistBank, species: Mapping[str, str], lineage: Mapping[str, Mapping[str, str]],
                    kingdom: str | None, ott_ncbi: Mapping[str, int | None] | None = None,
                    workers: int = 4) -> list[dict]:
    """Rows of col_species.tsv: each species name (to its genus) resolved in COL and crosswalked to NCBI.

    ``ott_ncbi`` gives the NCBI id OTT has for a species name, for comparison.
    """
    clb.info
    names = sorted(species)

    def one(name: str) -> dict:
        genus = species[name]
        ctx = {r: (lineage.get(genus) or {}).get(f"records_{r}", "") for r in LINEAGE_RANKS}
        row = resolve_name(clb, name, "species", ctx, kingdom, ott=False)
        ott_id = (ott_ncbi or {}).get(name)
        return {**{k: row.get(k, "") for k in SPECIES_COLUMNS}, "genus": genus, "species": name,
                "ott_ncbi_id": "" if ott_id is None else str(ott_id),
                "col_ncbi_vs_ott": compare_ids(row["col_ncbi_id"], "", "" if ott_id is None else str(ott_id))}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(one, names))
