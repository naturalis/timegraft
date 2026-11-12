"""Genus-level topology from static Open Tree of Life downloads.

Uses the OTT taxonomy (``taxonomy.tsv``) and the labelled synthetic tree
(``labelled_supertree.tre``) from files.opentreeoflife.org, instead of the
web API. Both are read once; the synthetic tree (about 2.5 million tips) is
held as flat parent and label arrays rather than as tree objects.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import numpy as np

# Taxa with these flags are passed over when an unflagged candidate exists.
DISFAVOURED_FLAGS = ("extinct", "barren", "hidden", "not_otu", "environmental", "incertae_sedis_inherited")


def _rows(path: str | Path):
    with open(path, encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            parts = line.rstrip("\n").split("\t|\t")
            if len(parts) < 7:
                continue
            yield parts


def _ncbi(sourceinfo: str) -> int | None:
    for source in sourceinfo.split(","):
        if source.startswith("ncbi:"):
            try:
                return int(source[5:])
            except ValueError:
                return None
    return None


class Taxonomy:
    """The parts of OTT needed to match genera and find their species."""

    def __init__(self, path: str | Path, genera: set[str], context: str | None = None):
        uids, parents, genus_uids = [], [], []
        self.ncbi_of: dict[int, int] = {}
        candidates: dict[str, list[tuple[int, str, int | None]]] = defaultdict(list)
        context_uid = None
        for uid, parent, name, rank, sourceinfo, _uniq, flags, *_ in _rows(path):
            uid_i = int(uid)
            uids.append(uid_i)
            parents.append(int(parent) if parent else -1)
            if rank in ("family", "superfamily", "order", "superorder", "subclass", "class", "subphylum", "phylum"):
                ncbi_id = _ncbi(sourceinfo)
                if ncbi_id:
                    self.ncbi_of[uid_i] = ncbi_id
            if rank == "genus":
                genus_uids.append(uid_i)
                if name in genera:
                    candidates[name].append((uid_i, flags, _ncbi(sourceinfo)))
            if context and name == context and context_uid is None:
                context_uid = uid_i
        if context and context_uid is None:
            raise ValueError(f"context taxon {context!r} not found in the taxonomy")
        self.parent = np.full(max(uids) + 1, -1, dtype=np.int64)
        self.parent[np.array(uids)] = np.array(parents)
        self.is_genus = np.zeros(len(self.parent), dtype=bool)
        self.is_genus[np.array(genus_uids, dtype=np.int64)] = True
        self.context_uid = context_uid

        self.genus: dict[str, dict] = {}
        for name in genera:
            options = [c for c in candidates.get(name, []) if self._within_context(c[0])]
            preferred = [c for c in options if not any(f in c[1] for f in DISFAVOURED_FLAGS)] or options
            if not preferred:
                self.genus[name] = {"ott_id": None, "ncbi_id": None, "status": "unmatched"}
            elif len(preferred) > 1:
                self.genus[name] = {"ott_id": None, "ncbi_id": None, "status": "ambiguous"}
            else:
                uid, _, ncbi = preferred[0]
                self.genus[name] = {"ott_id": uid, "ncbi_id": ncbi, "status": "matched"}

        # species (and lower) below each matched genus, with NCBI ids where known
        wanted = {m["ott_id"]: g for g, m in self.genus.items() if m["ott_id"] is not None}
        self.species: dict[str, list[tuple[int, int | None]]] = defaultdict(list)
        for uid, _parent, _name, rank, sourceinfo, _uniq, flags, *_ in _rows(path):
            if rank != "species" or "extinct" in flags:
                continue
            genus = self._ancestor_in(int(uid), wanted)
            if genus is not None:
                self.species[genus].append((int(uid), _ncbi(sourceinfo)))

    def _within_context(self, uid: int) -> bool:
        if self.context_uid is None:
            return True
        node, steps = uid, 0
        while node >= 0 and steps < 200:
            if node == self.context_uid:
                return True
            node = int(self.parent[node]) if node < len(self.parent) else -1
            steps += 1
        return False

    def _ancestor_in(self, uid: int, wanted: dict[int, str], max_steps: int = 6) -> str | None:
        node = int(self.parent[uid]) if uid < len(self.parent) else -1
        for _ in range(max_steps):
            if node < 0:
                return None
            if node in wanted:
                return wanted[node]
            node = int(self.parent[node])
        return None


def lookup_ids(taxonomy: "Taxonomy", genus: str, n_species: int = 12) -> list:
    """The genus's NCBI id (or None) followed by up to ``n_species`` species ids, spread over the list."""
    match = taxonomy.genus.get(genus) or {}
    species = sorted(n for _, n in taxonomy.species.get(genus, []) if n)
    if len(species) > n_species:
        step = len(species) / n_species
        species = [species[int(i * step)] for i in range(n_species)]
    return [match.get("ncbi_id") or 0] + species


def lineage(taxonomy: "Taxonomy", genus: str) -> list[tuple[int, int | None]]:
    """(OTT id, NCBI id or None) of every taxon enclosing a genus, nearest first."""
    match = taxonomy.genus.get(genus) or {}
    if not match.get("ott_id"):
        return []
    out, node = [], int(taxonomy.parent[match["ott_id"]])
    while node >= 0 and len(out) < 100:
        out.append((node, taxonomy.ncbi_of.get(node)))
        node = int(taxonomy.parent[node])
    return out


class SynthTree:
    """The labelled synthetic tree as flat arrays: parent index and label per node."""

    TOKEN = re.compile(r"[(),;]|[^(),;]+")

    def __init__(self, path: str | Path):
        text = Path(path).read_text().strip()
        parent: list[int] = []
        label: list[str] = []
        stack: list[int] = []
        current = -1
        last_closed = -1
        for match in self.TOKEN.finditer(text):
            token = match.group(0)
            if token == "(":
                node = len(parent)
                parent.append(stack[-1] if stack else -1)
                label.append("")
                stack.append(node)
                last_closed = -1
            elif token == ")":
                last_closed = stack.pop()
            elif token in ",;":
                last_closed = -1
            else:
                name = token.split(":")[0].strip().strip("'")
                if last_closed >= 0:
                    label[last_closed] = name
                    last_closed = -1
                else:
                    parent.append(stack[-1] if stack else -1)
                    label.append(name)
        self.parent = np.array(parent, dtype=np.int64)
        self.label = label
        self.index = {name: i for i, name in enumerate(label) if name}

    def _accumulate(self, values: list[int]) -> np.ndarray:
        # a parent is always read before its children, so one pass from the end suffices
        parent = self.parent.tolist()
        for i in range(len(parent) - 1, 0, -1):
            if parent[i] >= 0:
                values[parent[i]] += values[i]
        return np.array(values, dtype=np.int64)

    def tip_counts(self) -> np.ndarray:
        """Number of tips below every node (1 for a tip), in the order of ``parent``."""
        if getattr(self, "_tips", None) is None:
            self._tips = self._accumulate((~self.has_children()).astype(np.int64).tolist())
        return self._tips

    def has_children(self) -> np.ndarray:
        out = np.zeros(len(self.parent), dtype=bool)
        out[self.parent[self.parent >= 0]] = True
        return out

    def subtree_sizes(self) -> np.ndarray:
        """Number of nodes in the subtree of every node; nodes are numbered in preorder,
        so the subtree of ``n`` is ``range(n, n + size[n])``."""
        if getattr(self, "_sizes", None) is None:
            self._sizes = self._accumulate([1] * len(self.parent))
        return self._sizes

    def tip_genera(self, taxonomy: "Taxonomy", max_steps: int = 8) -> np.ndarray:
        """OTT id of the genus of every tip (itself if it is a genus), -1 where there is none."""
        out = np.full(len(self.parent), -1, dtype=np.int64)
        leaf = ~self.has_children()
        parent, is_genus = taxonomy.parent, taxonomy.is_genus
        for i in np.flatnonzero(leaf):
            name = self.label[i]
            if not name.startswith("ott"):
                continue
            uid = int(name[3:])
            for _ in range(max_steps):
                if uid < 0 or uid >= len(is_genus):
                    break
                if is_genus[uid]:
                    out[i] = uid
                    break
                uid = int(parent[uid])
        return out

    def path_to_root(self, node: int) -> list[int]:
        out = []
        while node >= 0:
            out.append(node)
            node = int(self.parent[node])
        return out

    def mrca(self, nodes: list[int]) -> int:
        paths = [self.path_to_root(n) for n in nodes]
        common = set(paths[0]).intersection(*map(set, paths[1:]))
        return next(n for n in paths[0] if n in common)


def genus_anchors(synth: SynthTree, taxonomy: Taxonomy, genera: list[str]) -> tuple[dict[str, int], dict[str, str]]:
    """The node of the synthetic tree that stands for each genus, and how it was found.

    See :func:`genus_topology_newick` for the rules and the status values.
    """
    anchor: dict[str, int] = {}
    status: dict[str, str] = {}
    for genus in genera:
        match = taxonomy.genus.get(genus, {"status": "unmatched"})
        if match["status"] != "matched":
            status[genus] = match["status"]
            continue
        node = synth.index.get(f"ott{match['ott_id']}")
        if node is not None:
            anchor[genus], status[genus] = node, "in_synthesis"
            continue
        present = [synth.index[f"ott{uid}"] for uid, _ in taxonomy.species.get(genus, []) if f"ott{uid}" in synth.index]
        if present:
            anchor[genus] = synth.mrca(present) if len(present) > 1 else int(synth.parent[present[0]])
            status[genus] = "attached_at_mrca"
            continue
        # last resort: the nearest taxonomic ancestor that is in the synthesis
        uid = int(taxonomy.parent[match["ott_id"]])
        while uid >= 0 and f"ott{uid}" not in synth.index:
            uid = int(taxonomy.parent[uid])
        if uid >= 0:
            anchor[genus], status[genus] = synth.index[f"ott{uid}"], f"attached_at_taxon:ott{uid}"
        else:
            status[genus] = "not_in_synthesis"
    return anchor, status


def genus_topology_newick(
    synth: SynthTree, taxonomy: Taxonomy, genera: list[str]
) -> tuple[str, dict[str, str]]:
    """Newick topology with one tip per genus, induced from the synthetic tree.

    A genus that is a node in the synthesis becomes a tip at that node. A
    genus that is not (not monophyletic, or only present through its
    species) is attached at the MRCA of its species in the synthesis, and
    failing that at its nearest taxonomic ancestor in the synthesis. Returns
    the Newick string and a status per genus: ``in_synthesis``,
    ``attached_at_mrca``, ``attached_at_taxon:<ott id>``, or why it is missing.
    """
    anchor, status = genus_anchors(synth, taxonomy, genera)
    if len(anchor) < 2:
        raise ValueError("fewer than two genera could be placed in the synthetic tree")

    children: dict[int, set[int]] = defaultdict(set)
    marked: set[int] = set()
    for node in anchor.values():
        for step in synth.path_to_root(node):
            if step in marked:
                break
            marked.add(step)
            parent = int(synth.parent[step])
            if parent >= 0:
                children[parent].add(step)
    for node in list(children):
        children[node] &= marked

    tips_at: dict[int, list[str]] = defaultdict(list)
    for genus, node in anchor.items():
        tips_at[node].append(genus)
    root = next(n for n in marked if int(synth.parent[n]) < 0 or int(synth.parent[n]) not in marked)

    def quote(name: str) -> str:
        return "'" + name.replace("'", "''") + "'"

    def build(node: int) -> str | None:
        parts = [quote(g) for g in sorted(tips_at.get(node, []))]
        sub = [s for s in (build(c) for c in sorted(children.get(node, ()))) if s]
        # a genus anchored here with nothing else below: it is the tip itself
        if not sub and len(parts) == 1:
            return parts[0]
        items = sub + parts
        if not items:
            return None
        if len(items) == 1:
            return items[0]  # collapse unary nodes
        return "(" + ",".join(items) + ")"

    import sys

    limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(limit, 100000))
    try:
        newick = build(root) + ";"
    finally:
        sys.setrecursionlimit(limit)
    return newick, status


def clade_richness(synth: SynthTree, anchor: dict[str, int], tree, taxonomy: "Taxonomy | None" = None) -> dict:
    """Size of the full clade in the synthetic tree for each internal node of ``tree``.

    ``tree`` is the genus topology built from the same ``anchor`` (genus to
    synthesis node), so each of its internal nodes is the MRCA of the anchors
    of its genera. With ``taxonomy`` the result is the number of genera in
    that clade (the genera of its tips), otherwise the number of tips, which
    are mostly species. Both include taxa that are not in the data, which is
    what a pure-birth model of the clade needs (#16).
    """
    if taxonomy is not None:
        genus_of = synth.tip_genera(taxonomy)
        sizes = synth.subtree_sizes()

        def size(m: int) -> int:
            values = genus_of[m : m + int(sizes[m])]
            return int(np.unique(values[values >= 0]).size)
    else:
        tips = synth.tip_counts()

        def size(m: int) -> int:
            return int(tips[m])

    paths: dict[int, list[int]] = {}

    def path(n: int) -> list[int]:
        if n not in paths:
            paths[n] = synth.path_to_root(n)
        return paths[n]

    def mrca(a: int, b: int) -> int:
        on_b = set(path(b))
        return next(n for n in path(a) if n in on_b)

    synth_node: dict = {}
    out: dict = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            synth_node[node] = anchor.get(node.taxon.label)
            continue
        found = [synth_node[c] for c in node.child_node_iter() if synth_node.get(c) is not None]
        if not found:
            synth_node[node] = None
            continue
        m = found[0]
        for other in found[1:]:
            m = mrca(m, other)
        synth_node[node] = m
        n = size(m)
        if n:
            out[node] = n
    return out


def crown_ages(taxonomy: Taxonomy, timetree, genera: list[str], max_species: int = 12) -> list[dict]:
    """Crown age per genus: the TimeTree MRCA of up to ``max_species`` of its species.

    This is a lower bound: the sampled species need not span the crown.
    Genera with fewer than two species known to TimeTree get no age.
    """
    picks_by_genus = {}
    for genus in genera:
        species = sorted(n for _, n in taxonomy.species.get(genus, []) if n)
        if len(species) < 2:
            continue
        step = max(1, len(species) // max_species)
        picks_by_genus[genus] = sorted(set(species[::step][:max_species]))
    if hasattr(timetree, "prefetch"):
        timetree.prefetch(picks_by_genus.values())
    rows = []
    for genus, picks in picks_by_genus.items():
        result = timetree._query(picks)
        if result["age"] is not None and len(result["found"]) >= 2:
            rows.append({"genus": genus, "crown_age": f"{result['age']:.4f}", "n_species": len(result["found"])})
    return rows
