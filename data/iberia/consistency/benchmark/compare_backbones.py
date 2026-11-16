"""Why the stem age of a genus differs between two backbones: the chain of causes per stem (#84).

    python compare_backbones.py A B OUT [--factor 1.10]

A and B are backbone folders (``backbone.tre``, ``node_ages.tsv`` and ``genera.tsv``, ``node_ages.tsv`` format 5
or later, which has ``age_source``, ``interpolated_from`` and ``query_genera``). For every genus whose stem age
(the age of the parent of its leaf) differs by more than ``--factor`` between A and B, the stem node of B is
followed back to what changed, node by node:

- ``topology``: the stem clade has other genera in A (a genus anchored or re-anchored elsewhere);
- for a node with its own calibration in both runs: ``identifier`` (a genus of the query was sent with another
  NCBI id, ``join_ncbi_id`` in ``genera.tsv``), ``higher_taxa`` (a query with enclosing taxa changed), ``picks``
  (other genera or species were picked for the query), ``timetree_answer`` (the same query, another answer),
  or ``fit`` (the same calibration, but the fit moved it: in one run it is kept, in the other fitted, or fitted
  differently);
- ``calibration_gained`` or ``calibration_lost``: the node is calibrated in one run only. This says what happened,
  not why, so a reason is added: ``no_answer`` (TimeTree has no age for any of the picks of the uncalibrated run;
  ``identifier`` instead if a genus of the node was sent with another id), ``identifier``, ``higher_taxa`` or
  ``picks`` (the two runs sent other queries), ``timetree_answer`` (the same query, another answer),
  ``unresolved`` (TimeTree's age is that of an ancestor, because its MRCA lies above the node), or ``fit`` (the
  other run had the calibration but its precedence rule dropped it as inconsistent with an ancestor; the fit
  keeps it). A stem whose only causes are a gained or lost calibration counts as not traced;
- for a node interpolated in both runs, what ``interpolate()`` places it by: its parent's age (followed if it
  changed, wherever the parent got it) and its floor, the oldest dated descendant in ``interpolated_from``
  (``anchors`` if it is another node, which is followed; followed too if only its age changed). If neither
  moved, ``richness`` if the clade richness of the node, its parent or its floor changed (the pure-birth
  fraction), or else ``topology`` (the number of undated nodes down to the floor differs);
- ``root``: the age of the root changed.

A step that is not a cause by itself (an interpolated node, or anchors that moved) leads to the nodes it
depends on. The chain ends at causes. Nodes are matched between the runs by their genera; ages within 0.1% count
as unchanged. Writes ``stem_causes.tsv`` (one row per changed stem: the causes, the chain, and whether it ends
in a cause) to OUT, and prints a summary.
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import dendropy

csv.field_size_limit(sys.maxsize)
SAME = 1.001  # ages within this factor count as unchanged
MAX_DEPTH = 12


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def factor(a, b):
    return max(a, b) / min(a, b) if a > 0 and b > 0 else float("inf")


class Backbone:
    def __init__(self, folder):
        folder = Path(folder)
        self.tree = dendropy.Tree.get(path=str(folder / "backbone.tre"), schema="newick",
                                      preserve_underscores=True, rooting="force-rooted")
        rows = {r["node"]: r for r in read(folder / "node_ages.tsv")}
        missing = {"age_source", "interpolated_from", "query_genera"} - set(next(iter(rows.values())))
        if missing:
            raise SystemExit(f"{folder}/node_ages.tsv lacks {sorted(missing)}: rebuild it with node_ages.tsv format 5")
        self.name, self.row, self.by_name, self.clade = {}, {}, {}, {}
        for i, node in enumerate(self.tree.preorder_node_iter()):
            name = f"n{i}"
            self.name[node], self.by_name[name] = name, node
            self.clade[node] = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
            if not node.is_leaf():
                self.row[node] = rows[name]
        self.by_clade = {c: n for n, c in self.clade.items()}
        self.leaf = {leaf.taxon.label: leaf for leaf in self.tree.leaf_node_iter()}
        self.genera = {r["genus"]: r for r in read(folder / "genera.tsv")}

    def age(self, node):
        return 0.0 if node.is_leaf() else float(self.row[node]["age"])

    def stem(self, genus):
        return self.leaf[genus].parent_node

    def describe(self, node):
        genera = sorted(self.clade[node])
        return f"{self.name[node]} ({len(genera)} genera: {';'.join(genera[:3])}{';...' if len(genera) > 3 else ''})"


def changed_identifiers(a, b, query_genera):
    """The genera of a query whose identifier sent to TimeTree (join_ncbi_id) differs between the runs."""
    genera = sorted({g for g in query_genera.split(";") if g and g in a.genera and g in b.genera
                     and a.genera[g].get("join_ncbi_id") != b.genera[g].get("join_ncbi_id")})
    return ", ".join(f"{g} {a.genera[g].get('join_ncbi_id') or '-'} to {b.genera[g].get('join_ncbi_id') or '-'}"
                     for g in genera[:5])


def calibrated(row):
    return row["age_source"] in ("timetree", "fit", "root") and bool(row["timetree_age"])


def explain(a, b, node_b, depth=0, seen=None):
    """(causes, chain) for why node_b of B has another age than its counterpart in A."""
    seen = set() if seen is None else seen
    if node_b in seen or depth > MAX_DEPTH:
        return [], []
    seen.add(node_b)
    clade = b.clade[node_b]
    node_a = a.by_clade.get(clade)
    here = b.describe(node_b)
    causes, chain = [], []
    if node_a is None:
        # another clade in A: say which genera moved, then go on against the node that stood in its place,
        # the most recent common ancestor in A of the same genera
        present = sorted(g for g in clade if g in a.leaf)
        if not present:
            return ["topology"], [f"{here}: none of its genera is in A"]
        mrca = a.tree.mrca(taxon_labels=present) if len(present) > 1 else a.leaf[present[0]].parent_node
        moved = []
        for g in sorted((clade ^ a.clade[mrca]) & set(a.genera) & set(b.genera)):
            ra_g, rb_g = a.genera[g], b.genera[g]
            if ra_g.get("anchor_node") != rb_g.get("anchor_node") or rb_g.get("placement_action") == "reanchored":
                moved.append(f"{g} ({'re-anchored' if rb_g.get('placement_action') == 'reanchored' else 'anchored'} "
                             f"{ra_g.get('anchor_node') or '-'} to {rb_g.get('anchor_node') or '-'})")
        detail = f"genera moved: {', '.join(moved[:5])}" if moved else "the clade's genera differ"
        causes.append("topology")
        chain.append(f"{here}: no node with these genera in A ({detail}); compared with {a.describe(mrca)}")
        node_a = mrca
        if node_a.is_leaf():
            return causes, chain
    more_causes, more_chain = _compare(a, b, node_a, node_b, here, depth, seen)
    return causes + more_causes, chain + more_chain


def _why_calibrated(a, b, with_cal, without, node_without, run_without, backbone_without):
    """(causes, notes) for why a node is calibrated in one run and not in the other.

    The build leaves a node uncalibrated where TimeTree has no age for any of its picks (query_ids empty),
    where TimeTree's age is that of an ancestor because its MRCA lies above the node (unresolved), or where
    the precedence rule drops the age as inconsistent with an ancestor (the fit keeps it). Behind the first
    two lies a query that differs between the runs: other identifiers, other picks or enclosing taxa.
    """
    causes, notes = [], []
    if not without["query_ids"]:
        # no query got an age; did the genera of the node send other identifiers?
        genera = ";".join(backbone_without.clade[node_without])
        changed = changed_identifiers(a, b, genera + ";" + with_cal["query_genera"])
        if changed:
            causes.append("identifier")
            notes.append(f"another identifier for its genera: {changed}")
        else:
            causes.append("no_answer")
            notes.append(f"TimeTree has no age for any of the picks of {run_without}")
        return causes, notes
    if without["query_ids"] != with_cal["query_ids"]:
        if "higher_taxa" in (without["timetree_attempt"], with_cal["timetree_attempt"]):
            causes.append("higher_taxa")
            notes.append("queried with enclosing taxa in one run")
        else:
            changed = changed_identifiers(a, b, without["query_genera"] + ";" + with_cal["query_genera"])
            if changed:
                causes.append("identifier")
                notes.append(f"its query sends another identifier: {changed}")
            else:
                causes.append("picks")
                notes.append(f"other picks for its query ({without['query_genera'] or '-'} in {run_without}, "
                             f"{with_cal['query_genera'] or '-'} in the other)")
    elif without["timetree_age"] != with_cal["timetree_age"]:
        causes.append("timetree_answer")
        notes.append(f"the same query, TimeTree {without['timetree_age'] or '-'} in {run_without}, "
                     f"{with_cal['timetree_age'] or '-'} in the other")
    if "unresolved" in without["source"] or "same_as_ancestor" in without["source"]:
        causes.append("unresolved")
        notes.append(f"TimeTree's age in {run_without} ({without['timetree_age']}) is that of an ancestor: "
                     f"its MRCA lies above the node")
    if "drop" in without["source"] or (without["timetree_age"] and "dropped" in without["precedence_source"]):
        causes.append("fit")
        notes.append(f"{run_without} has the calibration ({without['timetree_age']}) but dropped it as "
                     f"inconsistent with an ancestor ({without['source']})")
    return causes, notes


def _compare(a, b, node_a, node_b, here, depth, seen):
    """The causes of the age difference between node_a of A and node_b of B, as explained in the docstring."""
    clade = b.clade[node_b]
    ra, rb = a.row.get(node_a), b.row.get(node_b)
    if ra is None or rb is None:
        return [], []
    age_a, age_b = float(ra["age"]), float(rb["age"])
    if factor(age_a, age_b) <= SAME:
        return [], []
    step = f"{here}: {age_a:.2f} to {age_b:.2f} Ma"
    if node_b is b.tree.seed_node:
        return ["root"], [f"{step}, the root"]
    cal_a, cal_b = calibrated(ra), calibrated(rb)
    if cal_a and cal_b:
        if ra["query_ids"] != rb["query_ids"]:
            if rb["timetree_attempt"] == "higher_taxa" or ra["timetree_attempt"] == "higher_taxa":
                return ["higher_taxa"], [f"{step}; queried with other enclosing taxa"]
            changed = changed_identifiers(a, b, ra["query_genera"] + ";" + rb["query_genera"])
            if changed:
                return ["identifier"], [f"{step}; another identifier sent to TimeTree: {changed}"]
            return ["picks"], [f"{step}; other picks for its query ({ra['query_genera'] or '-'} to "
                               f"{rb['query_genera'] or '-'})"]
        if ra["timetree_age"] != rb["timetree_age"]:
            return ["timetree_answer"], [f"{step}; the same query, TimeTree {ra['timetree_age']} to {rb['timetree_age']}"]
        return ["fit"], [f"{step}; the same calibration ({rb['timetree_age']}), {ra['age_source']} in A, "
                         f"{rb['age_source']} in B"]
    if cal_a != cal_b:
        kind = "calibration_gained" if cal_b else "calibration_lost"
        causes, chain = [kind], [f"{step}; calibrated in {'B' if cal_b else 'A'} only (source {ra['source']} in A, "
                                 f"{rb['source']} in B)"]
        more, notes = _why_calibrated(a, b, rb if cal_b else ra, ra if cal_b else rb, node_a if cal_b else node_b,
                                      "A" if cal_b else "B", a if cal_b else b)
        return causes + more, chain + [f"  {n}" for n in notes]
    # interpolated in both. interpolate() places a node, from the root down, between the age of its parent
    # (itself perhaps interpolated) and that of its oldest dated descendant, its floor (or the tips), at a
    # fraction set by the richness of the node, its parent and its floor (or by the number of undated steps
    # down to the floor). So follow the floor and the parent, then compare those inputs.
    def floor_of(backbone, row):
        name = row["interpolated_from"].split(";")[1] if ";" in row["interpolated_from"] else ""
        return None if name in ("", "tips") else backbone.by_name.get(name)

    floor_a, floor_b = floor_of(a, ra), floor_of(b, rb)
    causes, chain = [], []
    same_floor = (floor_a is None and floor_b is None) or (
        floor_a is not None and floor_b is not None and a.clade[floor_a] == b.clade[floor_b])
    parent_b, parent_a = node_b.parent_node, node_a.parent_node
    parent_moved = factor(a.age(parent_a), b.age(parent_b)) > SAME
    floor_moved = same_floor and floor_b is not None and factor(a.age(floor_a), b.age(floor_b)) > SAME
    if not same_floor:
        chain.append(f"{step}; interpolated towards another dated descendant "
                     f"({ra['interpolated_from']} in A, {rb['interpolated_from']} in B)")
        causes.append("anchors")
        follow = [n for n in (floor_b,) if n is not None]
        if floor_b is None or a.by_clade.get(b.clade[floor_b]) is not None:
            follow += [n for n in (floor_a and b.by_clade.get(a.clade[floor_a]),) if n is not None]
    else:
        chain.append(f"{step}; interpolated between its parent and {rb['interpolated_from'].split(';')[1]}")
        follow = [floor_b] if floor_moved else []
    if parent_moved:
        follow.append(parent_b)
    for nxt in follow:
        c, s = explain(a, b, nxt, depth + 1, seen)
        causes += c
        chain += ["  " * (depth + 1) + x for x in s]
    if same_floor and not floor_moved and not parent_moved:
        inputs = [(node_a, node_b), (parent_a, parent_b)] + ([(floor_a, floor_b)] if floor_b is not None else [])
        changed = [(x, y) for x, y in inputs if x is not None and y is not None and not y.is_leaf()
                   and a.row[x]["clade_genera"] != b.row[y]["clade_genera"]]
        if changed:
            causes.append("richness")
            chain += [f"  clade richness of {b.describe(y)}: {a.row[x]['clade_genera']} to {b.row[y]['clade_genera']}"
                      for x, y in changed]
        else:
            causes.append("topology")
            chain.append("  the same parent age, floor and richness: the number of undated nodes down to the floor "
                         "differs")
    return causes, chain


ROOT_CAUSES = ("topology", "identifier", "higher_taxa", "picks", "timetree_answer", "unresolved", "no_answer", "fit",
               "calibration_gained", "calibration_lost", "richness", "root")
# a calibration gained or lost says what happened, not why: a stem is traced only with a cause besides these
NOT_A_REASON = {"calibration_gained", "calibration_lost"}


def traced(causes):
    """Whether a stem's causes explain it: a known cause besides a gained or lost calibration."""
    return bool(set(causes) & set(ROOT_CAUSES) - NOT_A_REASON)


def compare(a, b, threshold=1.10):
    """One row per genus whose stem age differs by more than threshold between backbones a and b."""
    rows = []
    for g in sorted(set(a.leaf) & set(b.leaf)):
        sa, sb = a.age(a.stem(g)), b.age(b.stem(g))
        if factor(sa, sb) <= threshold:
            continue
        causes, chain = explain(a, b, b.stem(g))
        found = sorted({c for c in causes if c in ROOT_CAUSES}, key=ROOT_CAUSES.index)
        rows.append({"genus": g, "stem_a": f"{sa:.2f}", "stem_b": f"{sb:.2f}", "factor": f"{factor(sa, sb):.3f}",
                     "causes": ";".join(found), "traced": "yes" if traced(found) else "no", "chain": " | ".join(chain)})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("out")
    ap.add_argument("--factor", type=float, default=1.10)
    args = ap.parse_args()
    rows = compare(Backbone(args.a), Backbone(args.b), args.factor)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "stem_causes.tsv", "w", newline="") as handle:
        w = csv.DictWriter(handle, list(rows[0]) if rows else ["genus"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n_traced = sum(1 for r in rows if r["traced"] == "yes")
    print(f"{len(rows)} stems differ by more than a factor {args.factor}; {n_traced} traced to a cause")
    print("causes:", dict(Counter(c for r in rows for c in r["causes"].split(";") if c)))
    for r in rows:
        if r["traced"] != "yes":
            print("untraced:", r["genus"], r["chain"][:300])


if __name__ == "__main__":
    main()
