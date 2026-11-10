"""Building one dated subtree per genus."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import dendropy

from . import distances
from .bionj import bionj
from .dating import crown_age, lsd2_ultrametric, mpl_ultrametric
from .rooting import root_genus
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
) -> tuple[dendropy.Tree, dict]:
    """Return the dated ingroup subtree and a provenance record.

    ``rooting`` holds ``rule`` (tolerant or strict) and ``fallback`` (mad or
    midpoint); the defaults are tolerant and mad.
    """
    rooting = rooting or {}
    if len(ingroup) < 2:
        raise ValueError("build_genus needs at least two units")
    method = marker_cfg["distance"]["method"]
    options = {k: v for k, v in marker_cfg["distance"].items() if k != "method"}
    items = prepare_sequences(ingroup + outgroup, method, marker_cfg.get("align", "none"))
    labels = [name for name, _ in items]
    dist = distances.matrix([seq for _, seq in items], method=method, **options)

    tree = bionj(labels, dist)
    tree, how, rooted = root_genus(tree, [n for n, _ in ingroup], [n for n, _ in outgroup],
                                   rule=rooting.get("rule", "tolerant"), fallback=rooting.get("fallback", "mad"))

    crown_cfg = marker_cfg["crown_age"]
    crown = crown_age(
        genus,
        stem_age,
        tree,
        known_crowns,
        method=crown_cfg["method"],
        rate=crown_cfg.get("rate"),
        fraction=crown_cfg.get("fraction"),
        cap_fraction=crown_cfg.get("cap_fraction", 0.9),
        suggested=timetree_crowns,
        max_ratio=crown_cfg.get("timetree_max_ratio", 3.0),
    )
    dating = marker_cfg.get("dating", {"method": "mpl"})
    if dating.get("method", "mpl") == "lsd2":
        tree = lsd2_ultrametric(
            tree, crown.age, seq_length=int(dating.get("seq_length", 650)),
            executable=dating.get("executable", "lsd2"),
        )
    else:
        mpl_ultrametric(tree, crown.age)
    provenance = {
        "genus": genus,
        "n_units": len(ingroup),
        "stem_age": f"{stem_age:.4f}",
        "crown_age": f"{crown.age:.4f}",
        "crown_source": crown.source,
        "timetree_crown": "" if crown.timetree is None else f"{crown.timetree:.4f}",
        "rooting": how,
        "violations": rooted["violations"],
        "misplaced": ";".join(rooted["misplaced"]),
        "dating": dating.get("method", "mpl"),
        "outgroup": ";".join(n for n, _ in outgroup),
        "status": "subtree",
    }
    return tree, provenance


def tree_to_newick(tree: dendropy.Tree) -> str:
    return tree.as_string(schema="newick", suppress_rooting=True, preserve_spaces=True).strip()
