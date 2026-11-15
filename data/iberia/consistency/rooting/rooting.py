"""Units on the wrong side of the root after the screening fix of #55 (#75).

    PYTHONPATH=src python data/iberia/consistency/rooting/rooting.py LEGACY DEFAULT OUT [--vsearch vsearch]

LEGACY and DEFAULT are the result folders (``results/iberia``) of two runs of the workflow through the genus
stage, on ``data/iberia/records.tsv.gz``: with the consistency layer switched off
(``benchmark/legacy.yaml``, as ``main`` ran before it) and at the default settings. Both need
``units/screened.tsv``, ``units/representatives.tsv``, ``genera/shard_*.tsv`` (the units of each genus and its
outgroup) and ``subtrees/shard_*.tsv`` (the provenance of each genus subtree).

A unit is on the wrong side of the root (counted in ``violations``) when it is a unit of the genus that the
root puts among the outgroup (``misplaced``), or a unit of the outgroup that it puts inside the genus. The
provenance names only the first kind, so for every genus whose wrong-side units may differ between the two
runs the unrooted genus tree is rebuilt as ``build_genus()`` builds it (the units aligned with MAFFT, as
``config/iberia.yaml`` sets, k2p distances, BIONJ; MAFFT must be on the path) and the minimum-violation edges of ``root_genus()`` are read from it; the edge chosen
is the one whose misplaced units are those recorded. The rebuilt counts must equal the recorded ones.

Writes to OUT:

- ``summary.tsv``: the counts of ``benchmark/results/end_to_end.tsv`` for both runs, and the comparison;
- ``genera.tsv``: every genus whose wrong-side units differ between the runs: units in each run, wrong-side
  units by kind, what changed (units added by the screening fix, the outgroup), the likely cause, and the
  crown age in both runs with its dependence on the root;
- ``units.tsv``: every unit on the wrong side in one run only: its kind, its screening decision in both runs,
  and its best identity to the genus's other units and to the outgroup's (for a misplaced unit), or to the
  genus and to its own genus (for an outgroup unit inside the genus).

Identities are VSEARCH global identities (``--usearch_global``, every hit from 50%) between the
representative sequences without gaps, as in the screening.
"""

import argparse
import csv
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from timegraft import distances
from timegraft.bionj import bionj
from timegraft.build import prepare_sequences
from timegraft.rooting import minimum_violation_sides

csv.field_size_limit(sys.maxsize)
# a crown age "depends strongly on the root" when its range over root positions (crown_age_min to
# crown_age_max) spans this factor: the definition that gives the 20 and 19 genera of benchmark/REPORT.md
FACTOR = 1.25


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    columns = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as handle:
        out = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", restval="")
        out.writeheader()
        out.writerows(rows)


class Run:
    def __init__(self, folder):
        folder = Path(folder)
        self.prov = {r["genus"]: r for f in sorted((folder / "subtrees").glob("shard_*.tsv")) for r in read(f)}
        self.screened = {r["unit_id"]: r for r in read(folder / "units/screened.tsv")}
        self.groups = defaultdict(lambda: {"in": [], "out": []})
        for f in sorted((folder / "genera").glob("shard_*.tsv")):
            for r in read(f):
                self.groups[r["genus"]][r["role"]].append((r["unit_id"], r["sequence"]))
        self.side, self.rebuilt = {}, {}

    def subtree(self, g):
        r = self.prov.get(g)
        return r if r and r["status"] == "subtree" else None

    def violations(self, g):
        r = self.subtree(g)
        return int(r["violations"]) if r and r["violations"] else 0

    def misplaced(self, g):
        r = self.subtree(g)
        return frozenset(r["misplaced"].split(";")) if r and r["misplaced"] else frozenset()

    def rebuild(self, genera, workers=2):
        """Rebuild the unrooted trees of ``genera`` (those with violations), in parallel."""
        todo = [g for g in genera if self.violations(g) and g not in self.side]
        with ProcessPoolExecutor(workers) as pool:
            for g, result in zip(todo, pool.map(sides, [self.groups[g] for g in todo])):
                self.rebuilt[g] = result

    def wrong_side(self, g):
        """(misplaced ingroup units, outgroup units inside the genus), from the rebuilt unrooted tree."""
        if g in self.side:
            return self.side[g]
        if not self.violations(g):
            self.side[g] = (frozenset(), frozenset())
            return self.side[g]
        count, ties = self.rebuilt.get(g) or sides(self.groups[g])
        if count != self.violations(g):
            raise AssertionError(f"{g}: rebuilt {count} violations, recorded {self.violations(g)}")
        recorded = self.misplaced(g)
        options = [inside for misplaced, inside in ties if misplaced == recorded]
        if not options:
            raise AssertionError(f"{g}: no minimum-violation edge has the recorded misplaced units")
        # the same set unless tied edges with the same misplaced units differ in the outgroup units inside
        self.side[g] = (recorded, frozenset().union(*options))
        return self.side[g]


def sides(group):
    """Minimum violations of a genus's unrooted tree, as ``build_genus()`` builds it, and for each tied edge
    the (misplaced ingroup units, outgroup units inside the genus)."""
    items = prepare_sequences(group["in"] + group["out"], "k2p", "mafft")
    dist = distances.matrix([s for _, s in items], method="k2p")
    tree = bionj([n for n, _ in items], dist)
    ingroup, outgroup = frozenset(n for n, _ in group["in"]), frozenset(n for n, _ in group["out"])
    count, ties = minimum_violation_sides(tree, ingroup, outgroup)
    return count, [(ingroup - side, outgroup & side) for _, side, _ in ties]


def strongly(row):
    lo, hi = (row or {}).get("crown_age_min", ""), (row or {}).get("crown_age_max", "")
    return bool(lo and hi and float(lo) > 0 and float(hi) / float(lo) >= FACTOR)


def identities(vsearch, sequences, queries, targets):
    """{(query, target): identity} for every hit of at least 50%."""
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        q, d, hits = Path(tmp) / "q.fa", Path(tmp) / "db.fa", Path(tmp) / "hits.tsv"
        q.write_text("".join(f">{u}\n{sequences[u].replace('-', '')}\n" for u in sorted(queries)))
        d.write_text("".join(f">{u}\n{sequences[u].replace('-', '')}\n" for u in sorted(targets)))
        subprocess.run([vsearch, "--usearch_global", str(q), "--db", str(d), "--id", "0.5", "--maxaccepts", "0",
                        "--maxrejects", "0", "--userfields", "query+target+id", "--userout", str(hits),
                        "--threads", "2", "--quiet"], check=True)
        for line in hits.read_text().splitlines():
            a, b, ident = line.split("\t")
            if a != b:
                out[(a, b)] = float(ident) / 100
    return out


def other_side(row):
    """Whether a wrong-side unit is closer to the side the root put it on than to its own: a misplaced unit
    to the outgroup, an outgroup unit to the genus."""
    return row["closer_to"] == ("outgroup" if row["kind"] == "misplaced" else "genus")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("legacy")
    ap.add_argument("default")
    ap.add_argument("out")
    ap.add_argument("--vsearch", default="vsearch")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    L, D = Run(args.legacy), Run(args.default)
    runs = {"legacy": L, "default": D}
    sequences, genus_of = {}, {}
    for folder in (args.legacy, args.default):
        for r in read(Path(folder) / "units/representatives.tsv"):
            sequences[r["unit_id"]], genus_of[r["unit_id"]] = r["sequence"], r["genus"]

    summary = []
    for name, run in runs.items():
        subtrees = [g for g in run.prov if run.subtree(g)]
        summary.append({
            "run": name,
            "genera with a subtree": len(subtrees),
            "units set aside by screening": sum(1 for r in run.screened.values()
                                                if r["decision"] in ("other_genus", "isolated")),
            "genera with units on the wrong side of the root": sum(1 for g in subtrees if run.violations(g)),
            "units on the wrong side of the root": sum(run.violations(g) for g in subtrees),
            "units of a genus among its outgroup": sum(len(run.misplaced(g)) for g in subtrees),
            f"crown ages whose range over root positions spans a factor {FACTOR}":
                sum(1 for g in subtrees if strongly(run.subtree(g))),
        })

    # genera whose wrong-side units differ in number or in the misplaced units of the genus. A genus with the
    # same count and the same misplaced units in both runs is not compared, even if the outgroup units inside
    # it could differ: that leaves the counts of end_to_end.tsv unchanged
    candidates = sorted(g for g in set(L.prov) | set(D.prov)
                        if (L.violations(g) or D.violations(g))
                        and (L.violations(g) != D.violations(g) or L.misplaced(g) != D.misplaced(g)))
    L.rebuild(candidates)
    D.rebuild(candidates)
    differ, only = [], {}
    for g in candidates:
        (lm, li), (dm, di) = L.wrong_side(g), D.wrong_side(g)
        if (lm, li) == (dm, di):
            continue
        differ.append(g)
        for name, (m, i), (om, oi) in (("default", (dm, di), (lm, li)), ("legacy", (lm, li), (dm, di))):
            only.update({(g, u): (name, "misplaced") for u in m - om})
            only.update({(g, u): (name, "outgroup_inside") for u in i - oi})

    # identities: a misplaced unit against the genus's other units and the outgroup; an outgroup unit inside
    # the genus against the genus's units and the other units of its own genus in that run
    compare = {}
    for (g, u), (name, kind) in only.items():
        run = runs[name]
        misplaced, inside = run.wrong_side(g)
        ingroup = {n for n, _ in run.groups[g]["in"]} - misplaced - {u}
        if kind == "outgroup_inside":
            compare[(g, u)] = (ingroup, {n for n, _ in run.groups[genus_of.get(u, "")]["in"]} - {u})
        else:
            compare[(g, u)] = (ingroup, {n for n, _ in run.groups[g]["out"]} - inside - {u})
    targets = set().union(*(a | b for a, b in compare.values())) if compare else set()
    queries = {u for _, u in only}
    ident = identities(args.vsearch, sequences, queries, targets | queries)

    unit_rows = []
    for (g, u), (name, kind) in sorted(only.items(), key=lambda x: (x[0][0], x[1][0], x[1][1], x[0][1])):
        a, b = compare[(g, u)]
        best_a = max(((ident.get((u, t), 0.0), t) for t in a), default=(0.0, ""))
        best_b = max(((ident.get((u, t), 0.0), t) for t in b), default=(0.0, ""))
        ld, dd = L.screened.get(u, {}), D.screened.get(u, {})
        row = {"genus": g, "unit_id": u, "wrong_side_in": name, "kind": kind, "unit_genus": genus_of.get(u, ""),
               "legacy_decision": ld.get("decision", ""), "default_decision": dd.get("decision", ""),
               "kept_by_the_fix": ld.get("decision") in ("other_genus", "isolated") and dd.get("decision") == "keep",
               "default_best_congener_identity": dd.get("best_congener_identity", ""),
               "default_best_other_genus": dd.get("best_other_genus", ""),
               "default_best_other_identity": dd.get("best_other_identity", ""),
               "identity_to_genus": f"{best_a[0]:.3f}" if best_a[1] else ""}
        if kind == "misplaced":
            row.update(identity_to_outgroup=f"{best_b[0]:.3f}" if best_b[1] else "",
                       nearest_outgroup_genus=genus_of.get(best_b[1], ""),
                       closer_to="outgroup" if best_b[0] > best_a[0] else "genus")
        else:
            row.update(identity_to_own_genus=f"{best_b[0]:.3f}" if best_b[1] else "",
                       closer_to="own_genus" if best_b[0] >= best_a[0] and best_b[1] else "genus")
        unit_rows.append(row)
    write(out / "units.tsv", unit_rows)

    by_genus = defaultdict(list)
    for r in unit_rows:
        by_genus[r["genus"]].append(r)
    genus_rows = []
    for g in differ:
        a, b = L.subtree(g) or {}, D.subtree(g) or {}
        units_a, units_b = {n for n, _ in L.groups[g]["in"]}, {n for n, _ in D.groups[g]["in"]}
        out_a, out_b = {n for n, _ in L.groups[g]["out"]}, {n for n, _ in D.groups[g]["out"]}
        new = [r for r in by_genus[g] if r["wrong_side_in"] == "default"]
        # units the fix kept, in this genus or in its outgroup, whether or not they are on the wrong side
        fix_units = sorted(u for u in (units_b - units_a) | (out_b - out_a)
                           if L.screened.get(u, {}).get("decision") in ("other_genus", "isolated")
                           and D.screened.get(u, {}).get("decision") == "keep")
        if not new:
            cause = "fewer wrong-side units in the default"
        elif any(r["kept_by_the_fix"] for r in new):
            cause = "a unit the screening fix kept is on the wrong side"
        elif fix_units:
            cause = "units the screening fix kept changed the tree"
        elif units_a != units_b or out_a != out_b:
            cause = "other changes to the genus or its outgroup"
        else:
            cause = "same units, other tree"  # not expected: the trees are deterministic
        crown_a, crown_b = a.get("crown_age", ""), b.get("crown_age", "")
        genus_rows.append({
            "genus": g,
            "wrong_side_only_in": "default" if D.violations(g) and not L.violations(g) else
                                  "legacy" if L.violations(g) and not D.violations(g) else "",
            "legacy_units": len(units_a), "default_units": len(units_b),
            "units_added": len(units_b - units_a), "units_removed": len(units_a - units_b),
            "outgroup_units_added": len(out_b - out_a), "outgroup_units_removed": len(out_a - out_b),
            "units_kept_by_the_fix": ";".join(fix_units),
            "legacy_violations": L.violations(g), "default_violations": D.violations(g),
            "legacy_misplaced": len(L.wrong_side(g)[0]), "legacy_outgroup_inside": len(L.wrong_side(g)[1]),
            "default_misplaced": len(D.wrong_side(g)[0]), "default_outgroup_inside": len(D.wrong_side(g)[1]),
            "new_wrong_side_units": ";".join(r["unit_id"] for r in new),
            "new_closer_to_the_other_side": sum(1 for r in new if other_side(r)),
            "cause": cause,
            "legacy_crown_age": crown_a, "default_crown_age": crown_b,
            "crown_factor": f"{max(float(crown_a), float(crown_b)) / min(float(crown_a), float(crown_b)):.3f}"
                            if crown_a and crown_b else "",
            "legacy_crown_source": a.get("crown_source", ""), "default_crown_source": b.get("crown_source", ""),
            "default_crown_age_min": b.get("crown_age_min", ""), "default_crown_age_max": b.get("crown_age_max", ""),
            "legacy_depends_strongly_on_root": strongly(a), "default_depends_strongly_on_root": strongly(b),
        })
    write(out / "genera.tsv", genus_rows)

    new_units = [r for r in unit_rows if r["wrong_side_in"] == "default"]
    factors = [float(r["crown_factor"]) for r in genus_rows if r["crown_factor"]]
    summary.append({
        "run": "comparison",
        "genera whose wrong-side units differ": len(genus_rows),
        "genera with wrong-side units in the default only": sum(1 for r in genus_rows if r["wrong_side_only_in"] == "default"),
        "genera with wrong-side units in legacy only": sum(1 for r in genus_rows if r["wrong_side_only_in"] == "legacy"),
        "units on the wrong side in the default only": len(new_units),
        "of which misplaced units of the genus": sum(1 for r in new_units if r["kind"] == "misplaced"),
        "of which outgroup units inside the genus": sum(1 for r in new_units if r["kind"] == "outgroup_inside"),
        "of which kept by the screening fix": sum(1 for r in new_units if r["kept_by_the_fix"]),
        "of which closer to the other side than to their own": sum(1 for r in new_units if other_side(r)),
        "units on the wrong side in legacy only": sum(1 for r in unit_rows if r["wrong_side_in"] == "legacy"),
        **{f"genera: {c}": n for c, n in Counter(r["cause"] for r in genus_rows).items()},
        "median crown factor, default against legacy": f"{statistics.median(factors):.3f}" if factors else "",
        "crown ages that differ by more than 10%": sum(1 for f in factors if f > 1.10),
        "crown ages that differ by more than a factor 1.5": sum(1 for f in factors if f > 1.5),
        "of these genera, default crowns that depend strongly on the root":
            sum(1 for r in genus_rows if r["default_depends_strongly_on_root"]),
        "of these genera, legacy crowns that depend strongly on the root":
            sum(1 for r in genus_rows if r["legacy_depends_strongly_on_root"]),
    })
    write(out / "summary.tsv", summary)
    for r in summary:
        print({k: v for k, v in r.items() if v != ""})


if __name__ == "__main__":
    main()
