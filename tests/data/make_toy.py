"""Generate the small synthetic data set used by the tests and the CI run.

Run from the repository root: python tests/data/make_toy.py
The output is committed, so this only needs rerunning if the design changes.
"""

import random
from pathlib import Path

OUT = Path(__file__).parent / "toy"
STOPS = {"TAA", "TAG"}  # invertebrate mitochondrial code (table 5)
BASES = "ACGT"

BACKBONE = "((Alpha:20,Beta:20):30,(Gamma:35,(Delta:15,Epsilon:15):20):15);\n"
CROWN_AGES = {"Alpha": 8.0}

# genus -> number of BINs; Epsilon has none and should be pruned
LAYOUT = {"Alpha": 4, "Beta": 3, "Gamma": 2, "Delta": 1}


def random_orf(rng, n_codons=104):
    codons = []
    while len(codons) < n_codons:
        codon = "".join(rng.choice(BASES) for _ in range(3))
        if codon not in STOPS:
            codons.append(codon)
    return "".join(codons)


def mutate(rng, seq, rate):
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    for i, codon in enumerate(codons):
        for pos in range(3):
            # third positions change ten times as often as first and second
            p = rate if pos == 2 else rate / 10
            if rng.random() < p:
                new = codon[:pos] + rng.choice([b for b in BASES if b != codon[pos]]) + codon[pos + 1 :]
                if new not in STOPS:
                    codon = new
        codons[i] = codon
    return "".join(codons)


def main():
    rng = random.Random(42)
    OUT.mkdir(parents=True, exist_ok=True)
    root = random_orf(rng)
    rows = []
    bin_number = 0
    for genus, n_bins in LAYOUT.items():
        ancestor = mutate(rng, root, 0.6)
        for _ in range(n_bins):
            bin_number += 1
            unit = f"BOLD:AAA{bin_number:04d}"
            bin_seq = mutate(rng, ancestor, 0.12)
            for r in range(3):
                rows.append((f"{unit.replace(':', '')}_{r}", unit, genus, mutate(rng, bin_seq, 0.003)))
    # a unit with no genus, to be reported for later placement
    rows.append(("BOLDAAZ0001_0", "BOLD:AAZ0001", "", mutate(rng, root, 0.5)))
    # a record with stop codons in every frame, to be rejected
    rows.append(("BAD_stop", "BOLD:AAA0005", "Beta", ("TAATAGATAAG" * 29)[:312]))
    # a record that is too short
    rows.append(("BAD_short", "BOLD:AAA0005", "Beta", root[:90]))

    with open(OUT / "records.tsv", "w") as handle:
        handle.write("record_id\tunit_id\tgenus\tsequence\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    (OUT / "backbone.tre").write_text(BACKBONE)
    with open(OUT / "crown_ages.tsv", "w") as handle:
        handle.write("genus\tcrown_age\n")
        for genus, age in CROWN_AGES.items():
            handle.write(f"{genus}\t{age}\n")


def substitute(rng, seq, n):
    seq = list(seq)
    for i in rng.sample(range(len(seq)), n):
        seq[i] = rng.choice([b for b in BASES if b != seq[i]])
    return "".join(seq)


def make_asvs():
    """Metabarcoding-style fragments of the reference records, plus a stranger."""
    rng = random.Random(7)
    records = {}
    with open(OUT / "records.tsv") as handle:
        next(handle)
        for line in handle:
            record_id, _, _, seq = line.rstrip("\n").split("\t")
            records[record_id] = seq
    start, end = 60, 270  # a 210 bp window inside the 312 bp barcode
    window = lambda rid: records[rid][start:end]
    asvs = {
        "asv_01": window("BOLDAAA0001_0"),                      # exact: matched
        "asv_02": substitute(rng, window("BOLDAAA0006_0"), 2),  # 99%: matched
        "asv_03": substitute(rng, window("BOLDAAA0008_0"), 13), # ~94%: placed near AAA0008
        "asv_04": "".join(rng.choice(BASES) for _ in range(210)),  # unrelated: unplaced
        "asv_05": substitute(rng, window("BOLDAAA0010_0"), 21), # 90%: placed near AAA0010
    }
    with open(OUT / "asvs.fasta", "w") as handle:
        for name, seq in asvs.items():
            handle.write(f">{name}\n{seq}\n")
    table = [
        ("estuary_1", "asv_01", 120), ("estuary_1", "asv_03", 30), ("estuary_1", "asv_04", 9),
        ("estuary_2", "asv_01", 60), ("estuary_2", "asv_02", 15), ("estuary_2", "asv_05", 22),
        ("marina_1", "asv_02", 200), ("marina_1", "asv_05", 41),
    ]
    with open(OUT / "asv_table.tsv", "w") as handle:
        handle.write("sample\tasv_id\treads\n")
        for row in table:
            handle.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
    make_asvs()
