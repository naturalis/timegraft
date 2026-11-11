"""Matching metabarcoding sequences (ASVs) to reference units, and placing the rest.

Search uses VSEARCH global alignment, which handles ASVs that cover only part
of the reference barcode (e.g. the Leray fragment of COI-5P).

Per ASV:
- ``matched``: best hit at or above ``match_identity``; the ASV counts as that unit.
- ``matched_ambiguous``: equally good best hits in several units; the ASV is
  attached at the most recent common ancestor of those units.
- ``placed``: best hit below ``match_identity`` but at or above
  ``place_identity``; the ASV is attached on the path from its nearest unit
  to the root, at the age implied by its distance to that unit.
- ``unplaced``: no hit at or above ``place_identity``.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import dendropy

from .seqio import write_fasta

PLACEMENT_COLUMNS = ("asv_id", "status", "tip", "unit_id", "genus", "identity", "attach_age")


def vsearch_hits(
    queries: list[tuple[str, str]],
    refs: list[tuple[str, str]],
    min_identity: float,
    query_cov: float = 0.9,
    maxaccepts: int = 20,
    threads: int = 1,
) -> dict[str, list[tuple[str, float]]]:
    """Hits per query as (reference id, identity in 0..1), best first."""
    if shutil.which("vsearch") is None:
        raise RuntimeError("vsearch not found on PATH")
    hits: dict[str, list[tuple[str, float]]] = defaultdict(list)
    if not queries or not refs:
        return hits
    with tempfile.TemporaryDirectory() as tmp:
        q, r, out = Path(tmp) / "q.fasta", Path(tmp) / "r.fasta", Path(tmp) / "hits.tsv"
        # neutral names: identifiers such as BOLD:AAA0001 are kept out of the FASTA headers
        write_fasta(q, [(f"q{i}", s.replace("-", "")) for i, (_, s) in enumerate(queries)])
        write_fasta(r, [(f"r{i}", s.replace("-", "")) for i, (_, s) in enumerate(refs)])
        subprocess.run(
            [
                "vsearch", "--usearch_global", str(q), "--db", str(r),
                "--id", str(min_identity), "--query_cov", str(query_cov),
                "--maxaccepts", str(maxaccepts), "--maxrejects", "64",
                "--userout", str(out), "--userfields", "query+target+id",
                "--threads", str(threads), "--quiet",
            ],
            check=True,
        )
        for line in out.read_text().splitlines():
            qi, ri, ident = line.split("\t")
            hits[queries[int(qi[1:])][0]].append((refs[int(ri[1:])][0], float(ident) / 100.0))
    for key in hits:
        hits[key].sort(key=lambda h: -h[1])
    return hits


def classify(
    hits: dict[str, list[tuple[str, float]]],
    asv_ids: list[str],
    record_unit: dict[str, str],
    unit_genus: dict[str, str],
    in_tree: set[str],
    match_identity: float,
    tie: float = 1e-6,
) -> list[dict]:
    """Decide per ASV whether it matches a unit, is to be placed, or neither.

    Only units that are tips of the tree count.
    """
    rows = []
    for asv in asv_ids:
        usable = [(record_unit[r], ident) for r, ident in hits.get(asv, []) if record_unit.get(r) in in_tree]
        if not usable:
            rows.append({"asv_id": asv, "status": "unplaced", "units": [], "identity": None})
            continue
        best = usable[0][1]
        units = sorted({u for u, ident in usable if best - ident <= tie})
        if best >= match_identity:
            status = "matched" if len(units) == 1 else "matched_ambiguous"
        else:
            status = "placed"
        genera = sorted({unit_genus.get(u, "") for u in units})
        rows.append(
            {"asv_id": asv, "status": status, "units": units, "identity": best,
             "genus": ";".join(g for g in genera if g)}
        )
    return rows


def _age(node: dendropy.Node) -> float:
    age = 0.0
    while not node.is_leaf():
        node = node.child_nodes()[0]
        age += node.edge.length or 0.0
    return age


def _mrca(tree: dendropy.Tree, labels: list[str]) -> dendropy.Node:
    nodes = [tree.find_node_with_taxon_label(label) for label in labels]
    paths = []
    for node in nodes:
        path = []
        while node is not None:
            path.append(node)
            node = node.parent_node
        paths.append(path)
    common = set(paths[0]).intersection(*map(set, paths[1:]))
    return next(n for n in paths[0] if n in common)


def attach(tree: dendropy.Tree, label: str, anchor: dendropy.Node, age: float) -> float:
    """Add tip ``label`` joining the path from ``anchor`` to the root at ``age``.

    If ``age`` is not above the anchor's own age, the tip joins the anchor
    directly (a polytomy). Returns the age at which it was attached.
    """
    namespace = tree.taxon_namespace
    tip = dendropy.Node(taxon=namespace.require_taxon(label=label))
    node, node_age = anchor, _age(anchor)
    if age <= node_age:
        node.add_child(tip)
        tip.edge.length = node_age
        return node_age
    while node.parent_node is not None:
        parent_age = node_age + (node.edge.length or 0.0)
        if age < parent_age:
            parent = node.parent_node
            joint = dendropy.Node()
            parent.remove_child(node)
            parent.add_child(joint)
            joint.edge.length = parent_age - age
            joint.add_child(node)
            node.edge.length = age - node_age
            joint.add_child(tip)
            tip.edge.length = age
            return age
        node, node_age = node.parent_node, parent_age
    # older than the root: join at the root
    node.add_child(tip)
    tip.edge.length = node_age
    return node_age


def place_all(
    tree: dendropy.Tree,
    rows: list[dict],
    height: str = "rate",
    rate: float | None = None,
    terminal_fraction: float = 0.5,
) -> list[dict]:
    """Attach placed and ambiguous ASVs to ``tree`` in place; return placement rows.

    ``height``: ``rate`` turns the distance to the nearest unit into an age
    with the per-lineage substitution rate (age = (1 - identity) / (2 * rate));
    ``terminal`` attaches at ``terminal_fraction`` of the nearest unit's
    terminal edge, for markers without a usable rate.
    """
    out = []
    for row in rows:
        record = {k: row.get(k) for k in ("asv_id", "status", "genus")}
        record["identity"] = f"{row['identity']:.4f}" if row["identity"] is not None else ""
        record["unit_id"] = ";".join(row["units"])
        record["tip"], record["attach_age"] = "", ""
        if row["status"] == "matched":
            record["tip"] = row["units"][0]
        elif row["status"] == "matched_ambiguous":
            anchor = _mrca(tree, row["units"])
            record["tip"] = row["asv_id"]
            record["attach_age"] = f"{attach(tree, row['asv_id'], anchor, 0.0):.4f}"
        elif row["status"] == "placed":
            anchor = _mrca(tree, row["units"]) if len(row["units"]) > 1 else tree.find_node_with_taxon_label(row["units"][0])
            if height == "rate":
                if not rate:
                    raise ValueError("placement height 'rate' needs a rate")
                age = (1.0 - row["identity"]) / (2.0 * rate)
            else:
                age = _age(anchor) + terminal_fraction * (anchor.edge.length or 0.0)
            record["tip"] = row["asv_id"]
            record["attach_age"] = f"{attach(tree, row['asv_id'], anchor, age):.4f}"
        out.append(record)
    return out


def attach_to_genus_stem(
    tree: dendropy.Tree,
    label: str,
    genus_units: list[str],
    fraction: float = 0.5,
) -> float | None:
    """Add tip ``label`` on the stem of a genus, at ``fraction`` of the way from its crown to its stem node.

    For units set aside as isolated within their genus (issue #38) that match
    nothing well enough to be placed elsewhere: they stay with their genus,
    as its earliest branch, but at an age set by the genus's crown and stem
    ages and not by their own distance. ``genus_units`` are the genus's units;
    those in the tree define the crown (a single one: its tip). Returns the
    attachment age, or None if no unit of the genus is in the tree.
    """
    labels = [u for u in genus_units if tree.find_node_with_taxon_label(u) is not None]
    if not labels:
        return None
    crown = _mrca(tree, labels) if len(labels) > 1 else tree.find_node_with_taxon_label(labels[0])
    if crown.parent_node is None:
        return None
    crown_age = _age(crown)
    stem_age = crown_age + (crown.edge.length or 0.0)
    return attach(tree, label, crown, crown_age + fraction * (stem_age - crown_age))

