"""Which node did TimeTree date? (#108, step 2)

The backbone is dated by asking TimeTree for the age of the most recent common ancestor of a few ids picked from
each child clade of a node. TimeTree answers from its own tree, which can differ from OToL's synthesis, so the
node it dates need not be the one meant.

TimeTree's tree refines the NCBI taxonomy: on the Iberian run, every named node TimeTree reported (2,158 of
2,164 whose lineages are known; the other six are NCBI clades added since TimeTree's NCBI release) is the NCBI
common ancestor of the ids it recognised. So the node it dated always lies within that common ancestor, and its
genera in the run can be compared with those of the backbone node (:func:`check`):

- ``match``: the same genera;
- ``narrower``: the common ancestor lacks genera of the backbone node, so TimeTree dated a node within it;
- ``wider``: the common ancestor also holds genera outside the backbone node. For a named answer, TimeTree's node
  is that common ancestor, so it dated a deeper node. For an unnamed one (a clade TimeTree made when resolving an
  NCBI polytomy) it may not have: :func:`probe` asks again with one such genus added. TimeTree reports its own id
  of the node it dated (``mrca_ttid``) for every answer; the same id means the genus is inside that node, so the
  answer dates a wider node (``wider``); another for every probe tried means it does not (``probed_outside``);
- ``unchecked``: a lineage is missing, or no probe could be made (TimeTree knows none of the candidate genera).

Genera without an NCBI id, or whose lineage is unknown, are left out of the comparison on both sides.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence

# the outcomes of check() and probe(); answers with these are not used as calibrations
REJECT = ("wider", "narrower")
OUTCOMES = ("match", "probed_outside", "wider", "narrower", "unchecked")


def common_ancestor(paths: Sequence[Sequence[str]]) -> str:
    """The deepest taxon shared by all ``paths`` (root first, each ending at the taxon itself), or ``""``."""
    prefix = ""
    for step in zip(*paths):
        if len(set(step)) > 1:
            break
        prefix = step[0]
    return prefix


def genus_lineages(genus_ids: Mapping[str, int | None],
                   lineage: Callable[[str], Sequence[str] | None]) -> dict[str, set[str]]:
    """Each genus's NCBI lineage as a set, for the genera whose id and lineage are known."""
    return {g: set(p) for g, i in genus_ids.items() if i and (p := lineage(str(i)))}


def check(found: Iterable[int], node_genera: Iterable[str], known: Mapping[str, set[str]],
          lineage: Callable[[str], Sequence[str] | None]) -> dict:
    """Compare the genera within the NCBI common ancestor of ``found`` with ``node_genera``.

    ``known`` gives the NCBI lineage of each genus of the run (:func:`genus_lineages`); ``lineage(taxon)`` the
    NCBI ids of a taxon from the root to the taxon itself, or None where unknown. Returns ``outcome``
    (``match``, ``narrower``, ``wider`` or ``unchecked``), ``ancestor`` (the NCBI common ancestor), and the
    genera ``outside`` it is wider by and ``missing`` from it.
    """
    paths = [lineage(str(i)) for i in found]
    if not paths or any(not p for p in paths):
        return {"outcome": "unchecked", "ancestor": "", "outside": set(), "missing": set()}
    ancestor = common_ancestor(paths)
    within = {g for g, path in known.items() if ancestor in path}
    node = set(node_genera) & set(known)
    outside, missing = within - node, node - within
    outcome = "wider" if outside else ("narrower" if missing else "match")
    return {"outcome": outcome, "ancestor": ancestor, "outside": outside, "missing": missing}


def probe(ids: Sequence[int], ttid, candidates: Sequence[str], genus_ids: Mapping[str, int],
          query: Callable[[list[int]], dict], probes: int = 3) -> tuple[str, list[str]]:
    """Whether TimeTree's node ``ttid`` for the query ``ids`` holds any of ``candidates`` (nearest first).

    ``query(ids)`` returns TimeTree's answer (with ``found`` and ``mrca_ttid``). Candidates TimeTree does not
    know are skipped; up to ``probes`` known ones are tried. Returns ``wider`` (a candidate is inside),
    ``probed_outside`` (every one tried is outside) or ``unchecked`` (none could be tried), and the genera tried.
    """
    if ttid in (None, "", 0):
        return "unchecked", []
    tried: list[str] = []
    for genus in candidates:
        extra = genus_ids.get(genus)
        if not extra or extra in ids:
            continue
        answer = query(sorted(set(ids) | {extra}))
        if extra not in (answer.get("found") or []):
            continue
        tried.append(genus)
        if answer.get("mrca_ttid") == ttid:
            return "wider", tried
        if len(tried) >= probes:
            break
    return ("probed_outside" if tried else "unchecked"), tried


def nearest_outside(node, allowed: set[str]) -> list[str]:
    """Genera in ``allowed`` outside the dendropy ``node``, nearest first: the sister clades, then those of the
    parent's sisters, and so on."""
    inside = {leaf.taxon.label for leaf in node.leaf_iter()}
    order: list[str] = []
    current = node
    while current.parent_node is not None:
        for sister in current.parent_node.child_node_iter():
            if sister is not current:
                order += sorted(leaf.taxon.label for leaf in sister.leaf_iter()
                                if leaf.taxon.label in allowed and leaf.taxon.label not in inside)
        current = current.parent_node
    return order
