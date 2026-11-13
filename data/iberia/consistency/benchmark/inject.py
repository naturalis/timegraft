"""Inject misplaced genera into a copy of the OTT taxonomy (#56).

    python inject.py TAXONOMY SYNTHESIS GENERA_TSV SEED OUT_TAXONOMY OUT_MANIFEST [N_PER_KIND] [merge|name|own]

The one real misplacement of the Iberian run, *Onchidella*, exists because OTT merged the genus with a
homonym elsewhere in the animals. This does the same to genera drawn at random: each is merged into a
genus-rank taxon that is a node of the synthesis and not in the run. That taxon takes the genus's name,
its NCBI ids (after its own) and its species as children; the genus's own taxon is left empty under
another name. Half the genera go to another phylum, half to another class of their own phylum.

Two variants are not part of the gate; they were added after its first results, to explain them.
With ``name``, only the name moves: the host keeps its own NCBI ids and children, so the merged-id
check has nothing to compare. With ``own``, the merge is as above, but the host's own species are
moved out of it (to the host's parent), so that the merged taxon holds the species of the genus only.
(The OTT entry of *Onchidella* holds the species of both, as the injection of the gate does.)

``GENERA_TSV`` is ``backbone/genera.tsv`` of the unperturbed default run: genera are drawn from those
with a place in the backbone and a single candidate in OTT. The pipeline itself is not touched.
"""

import csv
import random
import re
import sys

DISFAVOURED = ("extinct", "barren", "hidden", "not_otu", "environmental", "incertae_sedis", "unplaced",
               "major_rank_conflict", "merged", "inconsistent", "hybrid", "was_container")
SEP = "\t|\t"

taxonomy, synthesis, genera_tsv, seed, out_taxonomy, out_manifest = sys.argv[1:7]
n_per_kind = int(sys.argv[7]) if len(sys.argv) > 7 else 12
mode = sys.argv[8] if len(sys.argv) > 8 else "merge"
rng = random.Random(int(seed))

# 1. the taxonomy: parents, the ranks we need, and all genus rows
parent, rank_of, name_of, genus_rows = {}, {}, {}, {}
metazoa = None
with open(taxonomy) as handle:
    next(handle)
    for line in handle:
        uid, par, name, rank, sourceinfo, _uniq, flags, *_ = line.rstrip("\n").split(SEP)
        uid = int(uid)
        parent[uid] = int(par) if par else -1
        if rank in ("phylum", "class"):
            rank_of[uid], name_of[uid] = rank, name
        elif rank == "genus":
            genus_rows[uid] = (name, sourceinfo, flags)
        if name == "Metazoa" and metazoa is None:
            metazoa = uid


def lineage(uid):
    """(in Metazoa, phylum, class) of a taxon."""
    phylum = klass = ""
    inside = False
    while uid != -1:
        if uid == metazoa:
            inside = True
        if rank_of.get(uid) == "phylum" and not phylum:
            phylum = name_of[uid]
        if rank_of.get(uid) == "class" and not klass:
            klass = name_of[uid]
        uid = parent.get(uid, -1)
    return inside, phylum, klass


# 2. the run's genera, and the taxa they can be merged into
run = {r["genus"]: r for r in csv.DictReader(open(genera_tsv), delimiter="\t")}
eligible = sorted(g for g, r in run.items()
                  if r["ott_id"] and r["ott_candidates_in_context"] == "1" and g != "Onchidella"
                  and not r["status"].startswith(("unmatched", "ambiguous", "not_in_synthesis"))
                  and int(r["ott_id"]) in genus_rows)
nodes = {int(m) for m in re.findall(r"\)ott(\d+)", open(synthesis).read())}
hosts = {}  # (phylum, class) -> uids
for uid, (name, sourceinfo, flags) in genus_rows.items():
    if uid not in nodes or name in run or "ncbi:" not in sourceinfo or any(f in flags for f in DISFAVOURED):
        continue
    inside, phylum, klass = lineage(uid)
    if inside and phylum and klass:
        hosts.setdefault((phylum, klass), []).append(uid)
for uids in hosts.values():
    uids.sort()

# 3. draw genera and hosts
rng.shuffle(eligible)
chosen, used = [], set()
need = {"phylum": n_per_kind, "class": n_per_kind}
for genus in eligible:
    if not any(need.values()):
        break
    uid = int(run[genus]["ott_id"])
    _, phylum, klass = lineage(uid)
    if not phylum or not klass:
        continue
    kind = "phylum" if need["phylum"] >= need["class"] else "class"
    if kind == "phylum":
        options = [k for k in sorted(hosts) if k[0] != phylum]
    else:
        options = [k for k in sorted(hosts) if k[0] == phylum and k[1] != klass]
    if not options:
        continue
    host = rng.choice([u for u in hosts[rng.choice(options)] if u not in used])
    used.add(host)
    need[kind] -= 1
    _, host_phylum, host_class = lineage(host)
    chosen.append({"seed": seed, "genus": genus, "kind": kind, "ott_id": uid, "ott_phylum": phylum,
                   "ott_class": klass, "host_ott_id": host, "host_name": genus_rows[host][0],
                   "host_phylum": host_phylum, "host_class": host_class,
                   "records_class": run[genus]["records_class"], "records_phylum": run[genus]["records_phylum"]})

# 4. write the taxonomy with the merges
to_host = {c["ott_id"]: c["host_ott_id"] for c in chosen}
by_host = {c["host_ott_id"]: c for c in chosen}
with open(taxonomy) as src, open(out_taxonomy, "w") as dst:
    dst.write(next(src))
    for line in src:
        head = line.split(SEP, 2)
        uid = int(head[0])
        par = int(head[1]) if head[1] else -1
        if uid in to_host or uid in by_host or (par in to_host and mode != "name") or (par in by_host and mode == "own"):
            f = line.rstrip("\n").split(SEP)
            if uid in to_host:  # the genus's own taxon: left empty, under another name
                f[2] = f[2] + " zzmerged"
            elif uid in by_host:  # the host: takes the name and, after its own, the NCBI ids of the genus
                genus = by_host[uid]["genus"]
                own = [s for s in genus_rows[by_host[uid]["ott_id"]][1].split(",") if s.startswith("ncbi:")]
                sources = f[4].split(",")
                at = max(i for i, s in enumerate(sources) if s.startswith("ncbi:")) + 1
                f[2] = genus
                if mode != "name":
                    f[4] = ",".join(sources[:at] + own + sources[at:])
            elif par in to_host:  # a species of the genus: now a child of the host
                f[1] = str(to_host[par])
            else:  # a species of the host: out of the merged taxon, to the host's parent
                f[1] = str(parent[par])
            line = SEP.join(f) + "\n"
        dst.write(line)

with open(out_manifest, "w", newline="") as handle:
    writer = csv.DictWriter(handle, list(chosen[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(chosen)
print(f"seed {seed}: {len(chosen)} genera merged into hosts; {sum(len(v) for v in hosts.values())} hosts to choose from")
