"""Genus-level topology from the Open Tree of Life API (api.opentreeoflife.org, v3)."""

from __future__ import annotations

import json
import re
from typing import Callable, Iterable

import dendropy

from .webapi import ApiError, JsonCache, request_json

BASE = "https://api.opentreeoflife.org/v3"
BATCH = 250  # names per TNRS request


class OpenTree:
    def __init__(self, cache: JsonCache | None = None, fetch: Callable = request_json):
        self.cache = cache or JsonCache(None)
        self.fetch = fetch

    def match_genera(self, names: Iterable[str], context: str = "All life") -> dict[str, dict]:
        """Exact TNRS matches restricted to genus rank.

        Returns, per input name, a dict with ``ott_id``, ``ncbi_id`` (or None)
        and ``status``: ``matched``, ``unmatched`` or ``ambiguous``.
        """
        names = sorted(set(names))
        todo = [n for n in names if f"tnrs:{context}:{n}" not in self.cache]
        for start in range(0, len(todo), BATCH):
            chunk = todo[start : start + BATCH]
            result = self.fetch(
                f"{BASE}/tnrs/match_names",
                {"names": chunk, "context_name": context, "do_approximate_matching": False},
            )
            by_name = {r["name"]: r for r in result.get("results", [])}
            for name in chunk:
                self.cache.set(f"tnrs:{context}:{name}", _pick_genus(by_name.get(name)))
        return {n: self.cache.get(f"tnrs:{context}:{n}") for n in names}

    def induced_subtree(self, ott_ids: Iterable[int]) -> tuple[str, dict[int, str], dict[int, str]]:
        """Induced subtree of the synthetic tree, labelled with OTT ids.

        Returns the Newick string, the taxa left out because they are not in
        the synthetic tree (ott id -> reason), and the 'broken' taxa that are
        not monophyletic in the synthesis (ott id -> id of the node to attach
        them to).
        """
        ids = sorted({int(i) for i in ott_ids})
        missing: dict[int, str] = {}
        while True:
            try:
                result = self.fetch(
                    f"{BASE}/tree_of_life/induced_subtree", {"ott_ids": ids, "label_format": "id"}
                )
                break
            except ApiError as err:
                unknown = _unknown_ids(err)
                if err.status != 400 or not unknown:
                    raise
                for ott, reason in unknown.items():
                    missing[ott] = reason
                ids = [i for i in ids if i not in unknown]
                if len(ids) < 2:
                    raise RuntimeError("fewer than two genera are in the synthetic tree") from None
        broken = {
            int(k.replace("ott", "")): v for k, v in (result.get("broken") or {}).items() if k.startswith("ott")
        }
        return result["newick"], missing, broken


def _pick_genus(result: dict | None) -> dict:
    if not result:
        return {"ott_id": None, "ncbi_id": None, "status": "unmatched"}
    candidates = [
        m for m in result.get("matches", [])
        if m.get("taxon", {}).get("rank") == "genus" and not m.get("is_approximate_match")
    ]
    accepted = [m for m in candidates if not m.get("is_synonym")] or candidates
    if not accepted:
        return {"ott_id": None, "ncbi_id": None, "status": "unmatched"}
    ott_ids = {m["taxon"]["ott_id"] for m in accepted}
    if len(ott_ids) > 1:
        return {"ott_id": None, "ncbi_id": None, "status": "ambiguous"}
    taxon = accepted[0]["taxon"]
    ncbi = next((s.split(":", 1)[1] for s in taxon.get("tax_sources", []) if s.startswith("ncbi:")), None)
    return {"ott_id": int(taxon["ott_id"]), "ncbi_id": int(ncbi) if ncbi else None, "status": "matched"}


def _unknown_ids(err: ApiError) -> dict[int, str]:
    try:
        body = json.loads(err.body)
    except ValueError:
        return {}
    unknown = body.get("unknown") or {}
    out = {}
    for key, reason in unknown.items():
        match = re.fullmatch(r"ott(\d+)", key)
        if match:
            out[int(match.group(1))] = reason
    return out


def genus_topology(
    newick: str, ott_to_genus: dict[int, str], broken: dict[int, str]
) -> dendropy.Tree:
    """Relabel an OTT-labelled induced subtree with genus names.

    Broken genera are attached as an extra child of the node the API names
    for them, i.e. as part of a polytomy at their most recent common ancestor.
    """
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted", suppress_internal_node_taxa=True)
    by_label = {}
    for node in tree.preorder_node_iter():
        label = node.taxon.label if node.taxon is not None else node.label
        if label:
            by_label[label] = node
    namespace = tree.taxon_namespace
    for node in tree.leaf_nodes():
        ott = int(node.taxon.label.replace("ott", ""))
        node.taxon = namespace.require_taxon(label=ott_to_genus[ott])
    for ott, target in broken.items():
        parent = by_label.get(target)
        if parent is not None and ott in ott_to_genus:
            parent.add_child(dendropy.Node(taxon=namespace.require_taxon(label=ott_to_genus[ott])))
    for node in tree.preorder_node_iter():
        node.label = None
        if not node.is_leaf():
            node.taxon = None
    tree.suppress_unifurcations()
    keep = set(ott_to_genus.values())
    stray = [t for t in namespace if t.label not in keep]
    for taxon in stray:
        namespace.remove_taxon(taxon)
    return tree
