"""Building one dated subtree per genus."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import dendropy

from . import distances
from .bionj import bionj
from .dating import chronos_ultrametric, crown_age, lsd2_ultrametric, mean_root_to_tip, mpl_ultrametric
from .rooting import _root_split, copy_tree, root_genus, root_uncertainty
from .seqio import read_fasta, write_fasta


def align_mafft(items: list[tuple[str, str]], threads: int = 1) -> list[tuple[str, str]]:
    """Align with MAFFT; returns sequences in input order."""
    if shutil.which("mafft") is None:
        raise RuntimeError("mafft not found on PATH; install it or set align: none")
    with tempfile.TemporaryDirectory() as tmp:
        infile = Path(tmp) / "in.fasta"
        # neutral names so unit identifiers never upset the aligner
        write_fasta(infile, [(f"s{i}", seq.replace("-", "")) for i, (_, seq) in enumerate(items)])
        result = subprocess.run(
            ["mafft", "--auto", "--quiet", "--thread", str(threads), str(infile)],
            check=True,
            capture_output=True,
            text=True,
        )
        outfile = Path(tmp) / "out.fasta"
        outfile.write_text(result.stdout)
        aligned = dict(read_fasta(outfile))
    return [(name, aligned[f"s{i}"]) for i, (name, _) in enumerate(items)]


def prepare_sequences(items: list[tuple[str, str]], method: str, align: str) -> list[tuple[str, str]]:
    if not distances.NEEDS_ALIGNMENT.get(method, False):
        return items
    if align == "mafft":
        return align_mafft(items)
    if align == "none":
        lengths = {len(seq) for _, seq in items}
        if len(lengths) > 1:
            raise ValueError(
                f"distance method {method!r} needs aligned sequences but lengths differ; set align: mafft"
            )
        return items
    raise ValueError(f"unknown align option {align!r}")


def build_genus(
    genus: str,
    ingroup: list[tuple[str, str]],
    outgroup: list[tuple[str, str]],
    stem_age: float,
    marker_cfg: dict,
    known_crowns: dict[str, float],
    timetree_crowns: dict[str, float] | None = None,
    rooting: dict | None = None,
    rate: tuple[float | None, str] | None = None,
) -> tuple[dendropy.Tree, dict]:
    """Return the dated ingroup subtree and a provenance record.

    ``known_crowns`` holds, per genus, a point age or a mapping with ``age``,
    ``min`` and ``max`` (see ``dating.crown_age``). ``rate`` is the substitution
    rate for this genus and where it came from (see ``rates.Rates``); without
    it the marker's ``crown_age.rate`` is used.

    ``rooting`` holds ``rule`` (tolerant or strict), ``fallback`` (mad or
    midpoint) and ``uncertainty`` (whether to record how much the crown age
    and PD depend on the root, #34); the defaults are tolerant, mad and true.
    """
    rooting = rooting or {}
    if len(ingroup) < 2:
        raise ValueError("build_genus needs at least two units")
    method = marker_cfg["distance"]["method"]
    options = {k: v for k, v in marker_cfg["distance"].items() if k != "method"}
    items = prepare_sequences(ingroup + outgroup, method, marker_cfg.get("align", "none"))
    labels = [name for name, _ in items]
    dist = distances.matrix([seq for _, seq in items], method=method, **options)

    unrooted = bionj(labels, dist)
    in_names, out_names = [n for n, _ in ingroup], [n for n, _ in outgroup]
    uncertainty = rooting.get("uncertainty", True)
    before = copy_tree(unrooted) if uncertainty else None  # root_genus changes the tree it is given
    tree, how, rooted = root_genus(unrooted, in_names, out_names,
                                   rule=rooting.get("rule", "tolerant"), fallback=rooting.get("fallback", "mad"))

    crown_cfg = marker_cfg["crown_age"]
    genus_rate, rate_source = rate if rate is not None else (crown_cfg.get("rate"), "default")
    root_to_tip = mean_root_to_tip(tree)

    def choose_crown(t: dendropy.Tree):
        return crown_age(
            genus,
            stem_age,
            t,
            known_crowns,
            method=crown_cfg["method"],
            rate=genus_rate,
            fraction=crown_cfg.get("fraction"),
            cap_fraction=crown_cfg.get("cap_fraction", 0.9),
            suggested=timetree_crowns,
            timetree_bound=crown_cfg.get("timetree_bound", "max"),
            estimate=crown_cfg.get("estimate", "minimum"),
        )

    dating = marker_cfg.get("dating", {"method": "mpl"})

    if dating.get("method", "mpl") not in ("mpl", "lsd2", "chronos"):
        raise ValueError(f"unknown dating method {dating.get('method')!r}; use mpl, lsd2 or chronos")

    def date(t: dendropy.Tree, age: float) -> dendropy.Tree:
        if dating.get("method", "mpl") == "chronos":
            return chronos_ultrametric(
                t, age, model=dating.get("model", "correlated"), lam=float(dating.get("lambda", 1.0)),
                executable=dating.get("rscript", "Rscript"),
            )
        if dating.get("method", "mpl") == "lsd2":
            return lsd2_ultrametric(
                t, age, seq_length=int(dating.get("seq_length", 650)),
                executable=dating.get("executable", "lsd2"),
            )
        return mpl_ultrametric(t, age)

    by_rate_method = crown_cfg["method"] == "rate"
    crown = choose_crown(tree)
    chosen_split = _root_split(tree)
    tree = date(tree, crown.age)
    spread = {"root_ambiguity": "", "root_alternatives": "", "crown_age_min": "", "crown_age_max": "",
              "pd_rel_range": ""}
    if uncertainty:
        spread = _root_spread(before, in_names, out_names, rooted["misplaced"], tree, crown,
                              chosen_split, choose_crown, date)
    provenance = {
        "genus": genus,
        "n_units": len(ingroup),
        "stem_age": f"{stem_age:.4f}",
        "crown_age": f"{crown.age:.4f}",
        "crown_source": crown.source,
        "crown_lower": "" if crown.lower is None else f"{crown.lower:.4f}",
        "crown_upper": "" if crown.upper is None else f"{crown.upper:.4f}",
        "crown_upper_source": crown.upper_source,
        "barcode_age": "" if crown.barcode is None else f"{crown.barcode:.4f}",
        "timetree_crown": "" if crown.timetree is None else f"{crown.timetree:.4f}",
        "root_to_tip": f"{root_to_tip:.6f}",
        "rate": f"{genus_rate:.5f}" if by_rate_method and genus_rate else "",
        "rate_source": rate_source if by_rate_method and genus_rate else "",
        "rooting": how,
        "violations": rooted["violations"],
        "misplaced": ";".join(rooted["misplaced"]),
        **spread,
        "dating": dating.get("method", "mpl"),
        "outgroup": ";".join(n for n, _ in outgroup),
        "status": "subtree",
    }
    return tree, provenance


def tree_length(tree: dendropy.Tree) -> float:
    """Sum of branch lengths below the root (the PD of all tips)."""
    return sum(e.length or 0.0 for e in tree.postorder_edge_iter() if e.tail_node is not None)


def _root_spread(unrooted, ingroup, outgroup, misplaced, dated, crown, chosen_split, choose_crown, date) -> dict:
    """Crown ages and within-genus PD over the other plausible roots (issue #34).

    Crown age ranges are given only for crown ages taken from the barcodes
    (``rate`` sources, including those moved to a bound); table and fraction
    ages do not depend on the root. The PD range is relative to the PD of the chosen subtree and covers
    the alternatives with all units, not the one without the misplaced units.
    Alternatives whose crown age comes out at or below zero are skipped.
    """
    info = root_uncertainty(unrooted, ingroup, outgroup, misplaced)
    crowns, pds = [crown.age], [tree_length(dated)]
    done = {chosen_split}
    for label, alt in info["trees"]:
        full = label != "without_misplaced"
        if full:
            split = _root_split(alt)
            if split in done:
                continue
            done.add(split)
        alt_crown = choose_crown(alt)
        if alt_crown.age <= 0:
            # negative BioNJ branch lengths near this root: it cannot be dated, so it is not a plausible root
            continue
        crowns.append(alt_crown.age)
        if full:
            pds.append(tree_length(date(alt, alt_crown.age)))
    rai = info["root_ambiguity"]
    by_rate = crown.source.startswith(("rate", "midpoint"))
    return {
        "root_ambiguity": "" if math.isnan(rai) else f"{rai:.4f}",
        "root_alternatives": "" if info["root_alternatives"] is None else info["root_alternatives"],
        "crown_age_min": f"{min(crowns):.4f}" if by_rate else "",
        "crown_age_max": f"{max(crowns):.4f}" if by_rate else "",
        "pd_rel_range": f"{(max(pds) - min(pds)) / pds[0]:.4f}" if pds[0] > 0 else "",
    }


def tree_to_newick(tree: dendropy.Tree) -> str:
    return tree.as_string(schema="newick", suppress_rooting=True, preserve_spaces=True).strip()
