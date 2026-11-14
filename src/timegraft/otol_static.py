"""Genus-level topology from static Open Tree of Life downloads.

Uses the OTT taxonomy (``taxonomy.tsv``) and the labelled synthetic tree
(``labelled_supertree.tre``) from files.opentreeoflife.org, instead of the
web API. Both are read once; the synthetic tree (about 2.5 million tips) is
held as flat parent and label arrays rather than as tree objects.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from .timetree import EVIDENCE_COLUMNS

# Taxa with these flags are passed over when an unflagged candidate exists.
DISFAVOURED_FLAGS = ("extinct", "barren", "hidden", "not_otu", "environmental", "incertae_sedis_inherited")

# Ranks whose names are kept for the lineage of each genus (genera.tsv, #47).
LINEAGE_RANKS = ("family", "order", "class", "phylum")
# Taxa at or below these ranks are not kept by name (rank_name); everything above is, for lineages (#49).
LOW_RANKS = ("genus", "subgenus", "section", "subsection", "series", "species", "species group", "species subgroup",
             "subspecies", "variety", "varietas", "forma", "form", "subvariety", "no rank - terminal")


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


def _source_ids(sourceinfo: str, prefix: str) -> list[str]:
    """All ids of one source in an OTT ``sourceinfo`` (more than one where OTT merged records)."""
    return [s[len(prefix) + 1:] for s in sourceinfo.split(",") if s.startswith(prefix + ":")]


class Taxonomy:
    """The parts of OTT needed to match genera and find their species."""

    def __init__(self, path: str | Path, genera: set[str], context: str | None = None,
                 prefer: dict[str, int] | None = None):
        """``prefer`` maps a genus name to the OTT id to take when the name is otherwise ambiguous
        (for example from COL, #48); it is used only if that id is one of the candidates."""
        uids, parents, genus_uids = [], [], []
        self.ncbi_of: dict[int, int] = {}
        self.rank_name: dict[int, tuple[str, str]] = {}  # rank and name of taxa above genus
        candidates: dict[str, list[tuple[int, str, int | None]]] = defaultdict(list)
        sourceinfo_of: dict[int, str] = {}
        context_uid = None
        for uid, parent, name, rank, sourceinfo, _uniq, flags, *_ in _rows(path):
            uid_i = int(uid)
            uids.append(uid_i)
            parents.append(int(parent) if parent else -1)
            if rank in ("family", "superfamily", "order", "superorder", "subclass", "class", "subphylum", "phylum"):
                ncbi_id = _ncbi(sourceinfo)
                if ncbi_id:
                    self.ncbi_of[uid_i] = ncbi_id
            if rank not in LOW_RANKS:
                self.rank_name[uid_i] = (rank, name)
            if rank == "genus":
                genus_uids.append(uid_i)
                if name in genera:
                    candidates[name].append((uid_i, flags, _ncbi(sourceinfo)))
                    sourceinfo_of[uid_i] = sourceinfo
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
        # how each genus was matched, for genera.tsv; does not affect the match (#47)
        self.match_info: dict[str, dict] = {}
        for name in genera:
            options = [c for c in candidates.get(name, []) if self._within_context(c[0])]
            preferred = [c for c in options if not any(f in c[1] for f in DISFAVOURED_FLAGS)] or options
            if not preferred:
                self.genus[name] = {"ott_id": None, "ncbi_id": None, "status": "unmatched"}
            resolution = ""
            if len(preferred) > 1 and prefer and prefer.get(name) in {c[0] for c in preferred}:
                chosen = [c for c in preferred if c[0] == prefer[name]]
                resolution = "col"
            else:
                chosen = preferred
            if not chosen:
                self.genus[name] = {"ott_id": None, "ncbi_id": None, "status": "unmatched"}
            elif len(chosen) > 1:
                self.genus[name] = {"ott_id": None, "ncbi_id": None, "status": "ambiguous"}
            else:
                uid, _, ncbi = chosen[0]
                self.genus[name] = {"ott_id": uid, "ncbi_id": ncbi, "status": "matched"}
            info = {
                "ott_resolution": resolution,
                "ott_candidates": len(candidates.get(name, [])),
                "ott_candidates_in_context": len(options),
                "ott_candidate_ids": ";".join(str(c[0]) for c in preferred) if len(preferred) > 1 else "",
            }
            uid = self.genus[name]["ott_id"]
            if uid is not None:
                sourceinfo = sourceinfo_of.get(uid, "")
                info.update(
                    ott_flags=next(c[1] for c in preferred if c[0] == uid),
                    ott_sourceinfo=sourceinfo,
                    ott_ncbi_ids=";".join(_source_ids(sourceinfo, "ncbi")),
                    ott_worms_ids=";".join(_source_ids(sourceinfo, "worms")),
                    ott_gbif_ids=";".join(_source_ids(sourceinfo, "gbif")),
                )
            self.match_info[name] = info

        # species (and lower) below each matched genus, with NCBI ids where known
        wanted = {m["ott_id"]: g for g, m in self.genus.items() if m["ott_id"] is not None}
        self.species: dict[str, list[tuple[int, int | None]]] = defaultdict(list)
        self.species_ncbi: dict[str, int | None] = {}  # species name to its NCBI id, for comparison with COL
        self.species_uid: dict[str, int] = {}  # species name to its OTT id, for the species sampled in a run (#51)
        self.species_name: dict[int, str] = {}
        for uid, _parent, _name, rank, sourceinfo, _uniq, flags, *_ in _rows(path):
            if rank != "species" or "extinct" in flags:
                continue
            genus = self._ancestor_in(int(uid), wanted)
            if genus is not None:
                self.species[genus].append((int(uid), _ncbi(sourceinfo)))
                self.species_ncbi.setdefault(_name, _ncbi(sourceinfo))
                self.species_uid.setdefault(_name, int(uid))
                self.species_name[int(uid)] = _name

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


def lineage_names(taxonomy: "Taxonomy", genus: str) -> dict[str, str]:
    """Name of the nearest enclosing taxon at each of ``LINEAGE_RANKS`` in OTT, as ``ott_<rank>``."""
    out = {f"ott_{rank}": "" for rank in LINEAGE_RANKS}
    for uid, _ in lineage(taxonomy, genus):
        rank_name = taxonomy.rank_name.get(uid)
        if rank_name and rank_name[0] in LINEAGE_RANKS and not out[f"ott_{rank_name[0]}"]:
            out[f"ott_{rank_name[0]}"] = rank_name[1]
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


ANCHOR_RULES = ("majority", "mrca")
# genera.tsv columns describing how a genus that is not a node of the synthesis was anchored (#51)
ANCHOR_COLUMNS = ("anchor_rule", "anchor_fallback", "anchor_basis", "anchor_species", "anchor_species_inside",
                  "anchor_outside_share",
                  "anchor_outside_species", "anchor_mrca_node", "anchor_mrca_clade_genera")


def check_share(share: float) -> None:
    """The majority clade is only well defined for 0.5 <= share < 1: below a half, clades that are not
    nested can both pass; at 1 or more, none can."""
    if not 0.5 <= share < 1:
        raise ValueError(f"anchoring share must be at least 0.5 and less than 1, not {share}")


def majority_clade(synth: SynthTree, nodes: list[int], share: float = 0.5) -> tuple[int, int]:
    """The smallest clade of the synthetic tree holding more than ``share`` of ``nodes``, and how many it holds.

    Clades holding more than half of the same set are nested, so for
    ``share`` >= 0.5 the smallest is unique: the deepest node on the path
    from any of them to the root whose count passes the threshold.
    """
    check_share(share)
    counts: dict[int, int] = defaultdict(int)
    for node in set(nodes):
        for step in synth.path_to_root(node):
            counts[step] += 1
    n = len(set(nodes))
    passing = [node for node, c in counts.items() if c > share * n]
    sizes = synth.subtree_sizes()
    best = min(passing, key=lambda node: (int(sizes[node]), node))
    return best, counts[best]


def genus_anchors(synth: SynthTree, taxonomy: Taxonomy, genera: list[str], rule: str = "mrca",
                  sampled: dict[str, set[str]] | None = None, share: float = 0.5, min_species: int = 3,
                  min_sampled: int = 2, diagnostics: dict[str, dict[str, str]] | None = None,
                  ) -> tuple[dict[str, int], dict[str, str]]:
    """The node of the synthetic tree that stands for each genus, and how it was found.

    See :func:`genus_topology_newick` for the rules and the status values.

    A genus that is not a node of the synthesis is anchored on its species
    there. With ``rule`` ``mrca`` (the rule before #51) that is the MRCA of
    all of them, so a single stray species can drag the anchor, and the stem
    with it, deep into the tree. With ``rule`` ``majority`` (#51) it is the
    smallest clade holding more than ``share`` of them (``attached_at_majority``),
    where there are at least ``min_species``; with fewer, there is no
    majority to speak of and the MRCA is used. The species are those
    ``sampled`` in the run (names per genus) where at least ``min_sampled``
    of them are in the synthesis, otherwise all species of the genus in OTT.
    ``diagnostics``, if given, receives per genus the ``ANCHOR_COLUMNS``.
    """
    if rule not in ANCHOR_RULES:
        raise ValueError(f"anchoring rule must be one of {ANCHOR_RULES}, not {rule!r}")
    if rule == "majority":
        check_share(share)  # fail before any work, not at the first genus that needs it
    anchor: dict[str, int] = {}
    status: dict[str, str] = {}
    genus_of = None
    for genus in genera:
        match = taxonomy.genus.get(genus, {"status": "unmatched"})
        if match["status"] != "matched":
            status[genus] = match["status"]
            continue
        node = synth.index.get(f"ott{match['ott_id']}")
        if node is not None:
            anchor[genus], status[genus] = node, "in_synthesis"
            continue
        species = {uid for uid, _ in taxonomy.species.get(genus, []) if f"ott{uid}" in synth.index}
        basis = "all"
        if sampled and genus in sampled:
            mine = {taxonomy.species_uid[n] for n in sampled[genus] if n in taxonomy.species_uid} & species
            if len(mine) >= min_sampled:
                species, basis = mine, "sampled"
        present = sorted(synth.index[f"ott{uid}"] for uid in species)
        if present:
            mrca_node = synth.mrca(present) if len(present) > 1 else int(synth.parent[present[0]])
            anchor[genus], status[genus] = mrca_node, "attached_at_mrca"
            used, fallback = "mrca", ""
            if rule == "majority":
                if len(present) < min_species:
                    fallback = "too_few_species"
                else:
                    best, _ = majority_clade(synth, present, share)
                    if best != mrca_node:
                        anchor[genus], status[genus], used = best, "attached_at_majority", "majority"
                    else:
                        fallback = "no_smaller_majority"
            if diagnostics is not None:
                if genus_of is None:
                    genus_of = synth.tip_genera(taxonomy)
                sizes = synth.subtree_sizes()
                inside_nodes = set(range(anchor[genus], anchor[genus] + int(sizes[anchor[genus]])))
                outside = [uid for uid in sorted(species) if synth.index[f"ott{uid}"] not in inside_nodes]
                values = genus_of[mrca_node: mrca_node + int(sizes[mrca_node])]
                diagnostics[genus] = {
                    "anchor_rule": used,
                    "anchor_fallback": fallback,
                    "anchor_basis": basis,
                    "anchor_species": str(len(present)),
                    "anchor_species_inside": str(len(present) - len(outside)),
                    "anchor_outside_share": f"{len(outside) / len(present):.3f}",
                    "anchor_outside_species": ";".join(taxonomy.species_name.get(u, f"ott{u}") for u in outside[:10])
                    + (";..." if len(outside) > 10 else ""),
                    "anchor_mrca_node": synth.label[mrca_node],
                    "anchor_mrca_clade_genera": str(int(np.unique(values[values >= 0]).size)),
                }
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
    synth: SynthTree, taxonomy: Taxonomy, genera: list[str],
    anchors: tuple[dict[str, int], dict[str, str]] | None = None,
) -> tuple[str, dict[str, str]]:
    """Newick topology with one tip per genus, induced from the synthetic tree.

    A genus that is a node in the synthesis becomes a tip at that node. A
    genus that is not (not monophyletic, or only present through its
    species) is attached at the MRCA of its species in the synthesis, and
    failing that at its nearest taxonomic ancestor in the synthesis. Returns
    the Newick string and a status per genus: ``in_synthesis``,
    ``attached_at_mrca``, ``attached_at_majority`` (#51),
    ``attached_at_taxon:<ott id>``, or why it is missing.
    ``anchors`` gives the anchors and statuses to use instead, as returned by
    :func:`genus_anchors` and possibly corrected (#50).
    """
    anchor, status = anchors if anchors is not None else genus_anchors(synth, taxonomy, genera)
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

    out: dict = {}
    for node, m in synth_nodes(synth, anchor, tree).items():
        if node.is_leaf():
            continue
        n = size(m)
        if n:
            out[node] = n
    return out


def synth_nodes(synth: SynthTree, anchor: dict[str, int], tree) -> dict:
    """The node of the synthetic tree behind every node of ``tree`` (its index in ``synth``).

    A tip maps to the anchor of its genus, an internal node to the MRCA of
    the anchors of its genera. Nodes whose genera have no anchor are left out.
    """
    paths: dict[int, list[int]] = {}

    def path(n: int) -> list[int]:
        if n not in paths:
            paths[n] = synth.path_to_root(n)
        return paths[n]

    def mrca(a: int, b: int) -> int:
        on_b = set(path(b))
        return next(n for n in path(a) if n in on_b)

    out: dict = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            m = anchor.get(node.taxon.label)
        else:
            found = [out[c] for c in node.child_node_iter() if c in out]
            if not found:
                continue
            m = found[0]
            for other in found[1:]:
                m = mrca(m, other)
        if m is not None:
            out[node] = m
    return out


SUPPORT_KINDS = ("supported_by", "conflicts_with", "resolves", "partial_path_of", "terminal")


class _JSONStream:
    """Values from a JSON text, read incrementally from a file handle.

    Holds only a window of the text: each value is decoded with
    :meth:`json.JSONDecoder.raw_decode` once it is complete in the buffer, and
    the buffer is extended when it is not. Enough to walk the members of an
    object without loading the document; the members themselves are decoded
    whole.
    """

    CHUNK = 1 << 20
    _WS = re.compile(r"[ \t\n\r]*")

    def __init__(self, handle):
        self.handle = handle
        self.buf = ""
        self.pos = 0
        self.eof = False
        self.decoder = json.JSONDecoder()

    def _fill(self) -> bool:
        if self.eof:
            return False
        # grow the read with the pending text, so a long value is re-scanned a bounded number of times
        chunk = self.handle.read(max(self.CHUNK, len(self.buf) - self.pos))
        if not chunk:
            self.eof = True
            return False
        self.buf = self.buf[self.pos:] + chunk
        self.pos = 0
        return True

    def peek(self) -> str:
        """The next non-whitespace character, or ``""`` at the end."""
        while True:
            self.pos = self._WS.match(self.buf, self.pos).end()
            if self.pos < len(self.buf):
                return self.buf[self.pos]
            if not self._fill():
                return ""

    def expect(self, char: str) -> None:
        found = self.peek()
        if found != char:
            raise ValueError(f"malformed JSON: expected {char!r}, found {found!r}")
        self.pos += 1

    def value(self):
        self.peek()
        while True:
            try:
                obj, end = self.decoder.raw_decode(self.buf, self.pos)
            except json.JSONDecodeError:
                if self._fill():
                    continue
                raise
            # a number cut off at the end of the buffer decodes without error; make sure it is whole
            if end == len(self.buf) and self._fill():
                continue
            self.pos = end
            return obj

    def members(self):
        """(key, value) of the object that starts here, in file order."""
        self.expect("{")
        if self.peek() == "}":
            self.pos += 1
            return
        while True:
            key = self.value()
            self.expect(":")
            yield key, self
            if self.peek() == ",":
                self.pos += 1
            else:
                self.expect("}")
                return


def _top_level_member(handle, name: str):
    """(key, value) of the object under the top-level key ``name``, streamed (#72).

    Other top-level values are decoded and dropped one at a time; nothing is
    yielded when the key is absent or does not hold an object.
    """
    stream = _JSONStream(handle)
    for key, at in stream.members():
        if key == name and at.peek() == "{":
            for inner, value_at in at.members():
                yield inner, value_at.value()
        else:
            at.value()


def synthesis_support(path: str | Path, labels: set[str] | None = None) -> dict[str, dict[str, int]]:
    """Number of source trees per kind of annotation, per node of the synthetic tree.

    Read from the synthesis ``annotations.json``. The taxonomy is not a source
    tree there, so a node that is absent, or has no ``supported_by``, rests on
    the taxonomy alone. With ``labels`` only those nodes are kept.

    The file (61 MB, about 340,000 nodes for synthesis 16.1) is streamed rather
    than loaded, so only one node's annotation is held at a time (#72); as with
    :func:`json.load`, the last of duplicate keys wins.
    """
    out = {}
    with open(path, encoding="utf-8") as handle:
        for label, annotation in _top_level_member(handle, "nodes"):
            if labels is not None and label not in labels:
                continue
            if not isinstance(annotation, dict):
                out.pop(label, None)
                continue
            out[label] = {k: len(annotation[k]) for k in SUPPORT_KINDS if isinstance(annotation.get(k), dict)}
    return out


def genus_anchor_info(synth: SynthTree, taxonomy: Taxonomy, anchor: dict[str, int]) -> dict[str, dict[str, str]]:
    """Per genus: the synthesis node it is anchored at, and how much that clade holds (genera.tsv, #47).

    ``anchor_clade_genera`` is the number of genera in the clade of the anchor;
    for a genus attached at the MRCA of its species, a large value means one
    stray species pulled it deep (#51). ``n_species_ott`` and
    ``n_species_in_synthesis`` count its species in OTT and those of them that
    are tips or nodes of the synthetic tree.
    """
    genus_of = synth.tip_genera(taxonomy)
    sizes = synth.subtree_sizes()
    out = {}
    for genus, species in ((g, taxonomy.species.get(g, [])) for g in taxonomy.genus):
        row = {
            "anchor_node": "",
            "anchor_clade_genera": "",
            "n_species_ott": str(len(species)),
            "n_species_in_synthesis": str(sum(1 for uid, _ in species if f"ott{uid}" in synth.index)),
        }
        node = anchor.get(genus)
        if node is not None:
            values = genus_of[node : node + int(sizes[node])]
            row["anchor_node"] = synth.label[node]
            row["anchor_clade_genera"] = str(int(np.unique(values[values >= 0]).size))
        out[genus] = row
    return out


CROWN_COLUMNS = ("genus", "crown_age", "n_species", "query_ids", "crown_basis", "n_sampled_species",
                 "sampled_crown_age", "global_crown_age", *(f"timetree_{c}" for c in EVIDENCE_COLUMNS))
CROWN_BASES = ("sampled", "global")


def _spread(items: list, n: int) -> list:
    items = sorted(set(items))
    step = max(1, len(items) // n)
    return sorted(set(items[::step][:n]))


def crown_ages(taxonomy: Taxonomy, timetree, genera: list[str], max_species: int = 12,
               sampled: dict[str, set[str]] | None = None, basis: str = "sampled", min_sampled: int = 2,
               single_species: str = "none", species_ncbi: dict[str, int] | None = None) -> list[dict]:
    """Crown age per genus from TimeTree, and which species it is the crown of (``CROWN_COLUMNS``).

    The *global* crown is the TimeTree MRCA of up to ``max_species`` of the
    genus's species in OTT: the crown of the whole genus. In a regional run
    the genus subtree is built from the species sampled, and its root is the
    MRCA of those, which can be far younger (the holistic analysis: a
    *Ruditapes* of three units and 1.3% divergence against a global crown of
    69 Ma). With ``basis`` ``sampled`` (#54) the crown is therefore that of
    the species sampled in the run (``sampled``: species names per genus,
    found by name in OTT, or in ``species_ncbi``, e.g. from COL), where at
    least ``min_sampled`` of them have an NCBI id and TimeTree knows two;
    otherwise the global crown (``global_fallback``). A genus sampled with
    a single species has no TimeTree crown with ``single_species`` ``none``
    (``single_sampled_species``): its subtree is within one species, which
    the crown of the genus does not describe; ``global`` uses the global
    crown instead.

    Either age is a lower bound: the species queried need not span the
    crown. Genera with fewer than two species known to TimeTree get no age.
    Every genus with a decision gets a row; ``crown_age`` is empty where
    there is no TimeTree crown. Each row has both the sampled and the global
    age, the species queried and TimeTree's evidence for the chosen age.
    """
    if basis not in CROWN_BASES:
        raise ValueError(f"crown basis must be one of {CROWN_BASES}, not {basis!r}")
    if single_species not in ("none", "global"):
        raise ValueError(f"single_species must be 'none' or 'global', not {single_species!r}")
    names_to_ncbi = dict(taxonomy.species_ncbi)
    names_to_ncbi.update({k: v for k, v in (species_ncbi or {}).items() if v})
    global_picks, sampled_picks, n_sampled = {}, {}, {}
    for genus in genera:
        species = [n for _, n in taxonomy.species.get(genus, []) if n]
        if len(species) >= 2:
            global_picks[genus] = _spread(species, max_species)
        if sampled is not None and genus in sampled:
            ids = {names_to_ncbi[n] for n in sampled[genus] if names_to_ncbi.get(n)}
            n_sampled[genus] = len(ids)
            if len(ids) >= min_sampled:
                sampled_picks[genus] = _spread(list(ids), max_species)
    if hasattr(timetree, "prefetch"):
        timetree.prefetch(list(global_picks.values()) + (list(sampled_picks.values()) if basis == "sampled" else []))

    def age_of(picks):
        if not picks:
            return None
        result = timetree._query(picks)
        return result["age"] if result["age"] is not None and len(result["found"]) >= 2 else None

    rows = []
    for genus in genera:
        g_age = age_of(global_picks.get(genus))
        s_age = age_of(sampled_picks.get(genus)) if basis == "sampled" else None
        picks, age, how = None, None, "global"
        if basis == "sampled":
            if s_age is not None:
                picks, age, how = sampled_picks[genus], s_age, "sampled"
            elif n_sampled.get(genus) == 1 and single_species == "none":
                how = "single_sampled_species"
            elif g_age is not None:
                picks, age = global_picks[genus], g_age
                how = "global_fallback" if genus in n_sampled else "global"
        elif g_age is not None:
            picks, age = global_picks[genus], g_age
        if age is None and how != "single_sampled_species":
            continue
        found = timetree._query(picks)["found"] if picks else []
        row = {"genus": genus, "crown_age": "" if age is None else f"{age:.4f}", "n_species": len(found),
               "query_ids": "+".join(map(str, picks or [])), "crown_basis": how,
               "n_sampled_species": n_sampled.get(genus, ""),
               "sampled_crown_age": "" if s_age is None else f"{s_age:.4f}",
               "global_crown_age": "" if g_age is None else f"{g_age:.4f}"}
        if hasattr(timetree, "evidence"):
            row.update({f"timetree_{c}": v for c, v in
                        (timetree.evidence(picks) if picks else {c: "" for c in EVIDENCE_COLUMNS}).items()})
        rows.append(row)
    return rows
