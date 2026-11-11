"""Rates per taxon, their calibration, and the optional chronos dating."""

import shutil
import subprocess

import dendropy
import pytest

from timegraft.build import build_genus, tree_length
from timegraft.calibrate import calibrate, implied_rates, read_crown_table, trimmed_median
from timegraft.dating import chronos_ultrametric, is_ultrametric, mpl_ultrametric, root_age
from timegraft.rates import Rates, genus_taxonomy, read_rate_table


def test_rate_is_resolved_from_the_genus_upwards(tmp_path):
    table = tmp_path / "rates.tsv"
    table.write_text("rank\tname\trate\nall\tall\t0.004\nclass\tTeleostei\t0.0043\nphylum\tChordata\t0.005\n")
    assert read_rate_table(table) == {"class": {"Teleostei": 0.0043}, "phylum": {"Chordata": 0.005}}
    rates = Rates(0.0115, {"family": {"Gobiidae": 0.006}}, table)
    goby = {"genus": "Gobius", "family": "Gobiidae", "class": "Teleostei", "phylum": "Chordata"}
    assert rates.resolve(goby) == (0.006, "family:Gobiidae")
    assert rates.resolve({**goby, "family": "Sparidae"}) == (0.0043, "class:Teleostei")
    assert rates.resolve({"genus": "Ciona", "class": "Ascidiacea", "phylum": "Chordata"}) == (0.005, "phylum:Chordata")
    assert rates.resolve({"genus": "Nereis", "class": "Polychaeta"}) == (0.0115, "default")
    assert rates.resolve(None) == (0.0115, "default")
    assert Rates(None).resolve(goby) == (None, "default")  # markers without a rate
    with pytest.raises(ValueError):
        Rates(0.01, {"kingdom": {"Animalia": 0.01}})
    with pytest.raises(ValueError):
        Rates(0.01, {"class": {"Teleostei": 0.0}})


def test_genus_taxonomy_takes_the_most_common_value():
    records = [
        {"genus": "Gobius", "family": "Gobiidae", "class": "Teleostei"},
        {"genus": "Gobius", "family": "Gobiidae", "class": "Actinopterygii"},
        {"genus": "Gobius", "family": "", "class": "Teleostei"},
        {"genus": "", "family": "Gobiidae"},
    ]
    assert genus_taxonomy(records) == {
        "Gobius": {"genus": "Gobius", "family": "Gobiidae", "order": "", "class": "Teleostei", "phylum": ""}}


def test_trimmed_median_drops_the_tail_of_inflated_ages():
    good = [0.004, 0.0045, 0.005, 0.0055, 0.006]
    inflated = [0.0004, 0.0005, 0.0006]  # crown ages ten times too old
    centre, kept = trimmed_median(good + inflated)
    assert centre == pytest.approx(0.005) and sorted(kept) == good
    assert trimmed_median([0.01])[0] == 0.01


def test_calibration_per_class(tmp_path):
    def row(genus, distance, timetree="", stem=100.0, n=5):
        return {"genus": genus, "status": "subtree", "n_units": str(n), "stem_age": str(stem),
                "timetree_crown": str(timetree), "root_to_tip": str(distance)}

    provenance = [row(f"F{i}", 0.04 + 0.002 * i, 10.0) for i in range(10)]            # fishes, about 0.005
    provenance += [row(f"S{i}", 0.01 + 0.001 * i, 10.0) for i in range(10)]           # sharks, about 0.0015
    provenance += [row("Inflated", 0.05, 99.0), row("BeyondStem", 0.05, 120.0), row("NoAge", 0.05),
                   row("Single", 0.05, 10.0, n=1), row("Lit", 0.06, 500.0)]
    table = tmp_path / "crown.tsv"
    table.write_text("genus\tcrown_age\tcrown_age_min\nLit\t12\t\nOnlyMin\t\t145\n")
    crowns = read_crown_table(table)
    assert crowns == {"Lit": {"age": 12.0, "min": None, "max": None}, "OnlyMin": {"age": None, "min": 145.0, "max": None}}
    rates = implied_rates(provenance, {g: e["age"] for g, e in crowns.items() if e["age"]})
    assert "BeyondStem" not in rates and "NoAge" not in rates and "Single" not in rates
    assert rates["Lit"] == (pytest.approx(0.005), "table")  # the table age is used, not TimeTree's 500 Ma
    taxonomy = {f"F{i}": {"class": "Teleostei"} for i in range(10)}
    taxonomy.update({f"S{i}": {"class": "Elasmobranchii"} for i in range(10)})
    taxonomy.update({"Inflated": {"class": "Teleostei"}, "Lit": {"class": "Teleostei"}})
    rows = {(r["rank"], r["name"]): r for r in calibrate(rates, taxonomy, min_genera=10)}
    assert set(rows) == {("all", "all"), ("class", "Teleostei"), ("class", "Elasmobranchii")}
    fish = rows[("class", "Teleostei")]
    assert float(fish["rate"]) == pytest.approx(0.005, abs=3e-4) and fish["n_trimmed"] == 1  # Inflated dropped
    assert float(rows[("class", "Elasmobranchii")]["rate"]) == pytest.approx(0.00145, abs=1e-4)


def test_provenance_records_the_rate(genus_sequences):
    ingroup, outgroup = genus_sequences
    marker = {"align": "none", "distance": {"method": "k2p"},
              "crown_age": {"method": "rate", "rate": 0.0115, "cap_fraction": 0.9}, "dating": {"method": "mpl"}}
    _, default = build_genus("Genus", ingroup, outgroup, 400.0, marker, {}, rooting={"uncertainty": False})
    assert (default["rate"], default["rate_source"]) == ("0.01150", "default")
    assert float(default["crown_age"]) == pytest.approx(float(default["root_to_tip"]) / 0.0115, rel=1e-3)
    _, slow = build_genus("Genus", ingroup, outgroup, 400.0, marker, {}, rooting={"uncertainty": False},
                          rate=(0.00115, "class:Elasmobranchii"))
    assert (slow["rate"], slow["rate_source"]) == ("0.00115", "class:Elasmobranchii")
    assert float(slow["crown_age"]) == pytest.approx(10 * float(default["crown_age"]), rel=1e-3)
    assert slow["barcode_age"] == slow["crown_age"] and slow["crown_source"] == "rate"
    # TimeTree as an upper bound: the barcode age is kept in the record
    _, bound = build_genus("Genus", ingroup, outgroup, 400.0, marker, {}, {"Genus": 0.5},
                           rooting={"uncertainty": False})
    assert (bound["crown_source"], bound["crown_age"]) == ("rate_timetree_max", "0.5000")
    assert bound["barcode_age"] == default["crown_age"]


@pytest.fixture
def genus_sequences():
    import random

    rng = random.Random(3)

    def mutate(seq, n):
        seq = list(seq)
        for i in rng.sample(range(len(seq)), n):
            seq[i] = rng.choice([b for b in "ACGT" if b != seq[i]])
        return seq

    root = [rng.choice("ACGT") for _ in range(600)]
    clades = [mutate(root, 15), mutate(root, 15)]
    ingroup = [(f"G{i}", "".join(mutate(clades[i % 2], 5))) for i in range(6)]
    out = mutate(root, 60)
    outgroup = [(f"O{i}", "".join(mutate(out, 5))) for i in range(3)]
    return ingroup, outgroup


def _has_ape() -> bool:
    exe = shutil.which("Rscript")
    return bool(exe) and subprocess.run([exe, "-e", "library(ape)"], capture_output=True).returncode == 0


@pytest.mark.skipif(not _has_ape(), reason="needs R with the ape package")
def test_chronos_dates_to_the_given_root_age():
    newick = "(((A:0.02,B:0.03):0.02,(C:0.01,D:0.0):0.05):0.01,((E:0.04,F:0.02):0.01,G:0.06):0.02);"
    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted")
    dated = chronos_ultrametric(tree, 12.0)
    assert is_ultrametric(dated) and root_age(dated) == pytest.approx(12.0)
    assert sorted(l.taxon.label for l in dated.leaf_node_iter()) == list("ABCDEFG")
    # close to mean path length, which is why it is not the default
    mpl = mpl_ultrametric(dendropy.Tree.get(data=newick, schema="newick", rooting="force-rooted"), 12.0)
    assert tree_length(dated) == pytest.approx(tree_length(mpl), rel=0.25)
    with pytest.raises(ValueError):
        chronos_ultrametric(tree, 12.0, model="strict")


def test_chronos_falls_back_on_mean_path_length_for_small_trees():
    tree = dendropy.Tree.get(data="((A:0.1,B:0.3):0.2,C:0.25);", schema="newick", rooting="force-rooted")
    dated = chronos_ultrametric(tree, 12.0, executable="no-such-rscript")
    assert is_ultrametric(dated) and root_age(dated) == pytest.approx(12.0)
