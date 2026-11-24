"""Which node did TimeTree date? (#108, step 2)

The backbone is dated by asking TimeTree for the age of the most recent common ancestor of a few ids picked from
each child clade of a node. TimeTree answers from its own tree, which can differ from OToL's synthesis, so the
node it dates need not be the one meant.

TimeTree's tree refines the NCBI taxonomy: on the Iberian run, every named node TimeTree reported (2,158 of
2,164 whose lineages are known; the other six are NCBI clades added since TimeTree's NCBI release) is the NCBI
common ancestor of the ids it recognised. The genera of the run within a taxon (by their NCBI lineages) are then
compared with those of the backbone node (:func:`check`):

- a **named** answer is checked against the taxon TimeTree reports (``mrca_ncbi``), which is the node it dated;
- an **unnamed** answer (a clade TimeTree made when resolving an NCBI polytomy) lies within the NCBI common
  ancestor of the ids TimeTree recognised, and is checked against that.

The outcomes:

- ``narrower``: the taxon lacks genera of the backbone node, so TimeTree dated a node within it. This holds
  whether or not the taxon also holds genera outside the node;
- ``wider``: the taxon holds genera outside the node, and none of the node is missing. For a named answer that
  settles it; an unnamed one may still not hold them, so :func:`probe` asks again with one such genus added.
  TimeTree reports its own id of the node it dated (``mrca_ttid``) for every answer: the same id means the genus
  is inside that node (``wider``); another for every probe tried means it is outside (``probed_outside``);
- ``match``: the same genera;
- ``unchecked``: a lineage is missing, or no probe could be made.

An unnamed answer that passes (``match``, ``probed_outside``) can still date a node that leaves out genera of the
backbone node: the NCBI taxon holds them, but TimeTree's clade within it may not. :func:`probe` therefore also
asks with genera of the node that none of the recognised ids stand for; one outside TimeTree's node makes the
answer ``narrower``.

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
          lineage: Callable[[str], Sequence[str] | None], reported=None) -> dict:
    """Compare the genera of the run within TimeTree's node with ``node_genera``.

    ``known`` gives the NCBI lineage of each genus of the run (:func:`genus_lineages`); ``lineage(taxon)`` the
    NCBI ids of a taxon from the root to the taxon itself, or None where unknown. ``reported`` is the NCBI taxon
    TimeTree reports for the answer (``mrca_ncbi``), or 0 or None for an unnamed node. Returns ``outcome``
    (``match``, ``narrower``, ``wider`` or ``unchecked``), ``ancestor`` (the taxon compared: the reported one,
    or for an unnamed answer the NCBI common ancestor of ``found``), and the genera ``outside`` it is wider by
    and ``missing`` from it.
    """
    none = {"outcome": "unchecked", "ancestor": "", "outside": set(), "missing": set()}
    paths = [lineage(str(i)) for i in found]
    if not paths or any(not p for p in paths):
        return none
    if reported:
        ancestor = str(reported)
        # TimeTree's taxon must hold every recognised id; where today's NCBI disagrees (another NCBI release),
        # the answer cannot be checked safely
        if not lineage(ancestor) or any(ancestor not in p for p in paths):
            return {**none, "ancestor": ancestor}
    else:
        ancestor = common_ancestor(paths)
    within = {g for g, path in known.items() if ancestor in path}
    node = set(node_genera) & set(known)
    outside, missing = within - node, node - within
    outcome = "narrower" if missing else ("wider" if outside else "match")
    return {"outcome": outcome, "ancestor": ancestor, "outside": outside, "missing": missing}


def represented(found: Iterable[int], node_genera: Iterable[str], genus_ids: Mapping[str, int],
                lineage: Callable[[str], Sequence[str] | None]) -> set[str]:
    """The genera of the node that a recognised id stands for: the genus itself, or one of its species."""
    ids = {str(i) for f in found for i in (lineage(str(f)) or [str(f)])}
    return {g for g in node_genera if genus_ids.get(g) and str(genus_ids[g]) in ids}


def probe(ids: Sequence[int], ttid, candidates: Sequence[str], genus_ids: Mapping[str, int],
          query: Callable[[list[int]], dict], probes: int = 3, stop_when_inside: bool = True
          ) -> tuple[str, list[tuple[str, bool]]]:
    """Whether genera are inside TimeTree's node ``ttid`` for the query ``ids``, nearest ``candidates`` first.

    ``query(ids)`` returns TimeTree's answer (with ``found`` and ``mrca_ttid``). Candidates TimeTree does not
    know, or whose answer has no node id, are skipped; up to ``probes`` others are tried. With
    ``stop_when_inside`` (genera outside the backbone node) the first one inside ends the probe: ``wider``;
    otherwise (genera of the node not yet stood for) the first one outside does: ``narrower``. Returns that
    outcome, ``passed`` where every genus tried came out as it should, or ``unchecked`` where none could be
    tried, and the genera tried with whether each is inside.
    """
    if ttid in (None, "", 0):
        return "unchecked", []
    tried: list[tuple[str, bool]] = []
    for genus in candidates:
        extra = genus_ids.get(genus)
        if not extra or extra in ids:
            continue
        answer = query(sorted(set(ids) | {extra}))
        if extra not in (answer.get("found") or []) or answer.get("mrca_ttid") in (None, "", 0):
            continue
        inside = answer.get("mrca_ttid") == ttid
        tried.append((genus, inside))
        if inside and stop_when_inside:
            return "wider", tried
        if not inside and not stop_when_inside:
            return "narrower", tried
        if len(tried) >= probes:
            break
    return ("passed" if tried else "unchecked"), tried


def interleave(node) -> list[str]:
    """The genera of the dendropy ``node``, one from each child in turn, so that probes spread over the node."""
    lists = [sorted(leaf.taxon.label for leaf in child.leaf_iter()) for child in node.child_node_iter()]
    out: list[str] = []
    for k in range(max((len(x) for x in lists), default=0)):
        out += [x[k] for x in lists if k < len(x)]
    return out


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
