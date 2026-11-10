# How timegraft works

This page explains the method without code. For settings and commands, see the [technical documentation](../technical/configuration.md).

## The problem in one paragraph

Metabarcoding tells us which organisms are in a sample. To measure how much evolutionary history a sample holds, or how different two samples are, we need a tree that connects those organisms and whose branch lengths mean something, ideally time. DNA barcodes are short (about 650 bases for COI, a few hundred for ITS). They are good at telling close relatives apart, but they lose track of how far apart distant groups are: after enough time, mutations pile up on the same positions and the differences stop growing (saturation). Older approaches tried to build the whole tree from barcodes anyway, and failed at the deep end.

## The idea: let each source do what it is good at

```mermaid
flowchart LR
    subgraph external["Published knowledge"]
        otol["Open Tree of Life:<br/>how genera are related"]
        timetree["TimeTree:<br/>when lineages split"]
    end
    subgraph barcodes["Barcode libraries"]
        bold["BOLD (COI) or UNITE (ITS):<br/>reference sequences, grouped<br/>into BINs or species hypotheses"]
    end
    otol --> backbone["Dated backbone:<br/>one tip per genus,<br/>ages in millions of years"]
    timetree --> backbone
    bold --> genus["Within each genus:<br/>a small dated tree<br/>of its BINs"]
    backbone --> graft["Graft each genus tree<br/>onto its tip in the backbone"]
    genus --> graft
    graft --> tree["One dated tree with<br/>a BIN at every tip"]
```

- **Above genus level**, the tree comes from published knowledge: the branching order from the Open Tree of Life, which synthesises thousands of published phylogenies, and the ages of the branching points from TimeTree, which synthesises published divergence times.
- **Within each genus**, the barcodes decide: which BINs are most closely related, and how long ago they split. This is the range in which barcodes are reliable.
- **Grafting** joins the two: each genus tip in the backbone is replaced by the genus's own small tree.

Because everything is expressed in millions of years, the pieces fit together without any rescaling, and a tree for COI and a tree for ITS measure diversity in the same unit.

## Step 1: reference units

Barcode libraries hold many sequences per species. BOLD groups them into BINs, clusters that correspond roughly to species; UNITE does the same for fungi with species hypotheses. timegraft keeps one representative sequence per BIN (the one most similar to the others), after removing sequences that are too short, too ambiguous, or that cannot be read as a protein-coding gene. Each BIN is assigned to the genus that most of its records are identified as.

## Step 2: the dated backbone

The backbone is a tree of genera with known ages:

1. All genera in the data are looked up in the Open Tree of Life, which gives how they are related.
2. For each branching point, TimeTree is asked when the lineages on either side split. TimeTree does not know every genus, so it is also asked about species within the genera and, as a last resort, about the families or orders they belong to.
3. Branching points that TimeTree cannot date are given ages evenly spread between the nearest dated points above and below.

Each genus in the backbone then has a **stem age**: when it split from its closest relatives.

## Step 3: a small tree per genus

```mermaid
flowchart TB
    bins["The genus's BINs<br/>+ ten BINs from<br/>the closest related genera"] --> distances["How different is each pair?<br/>(barcode distances)"]
    distances --> join["Join the most similar<br/>step by step into a tree"]
    join --> root["Find the root:<br/>where the related genera<br/>attach"]
    root --> drop["Remove the related genera"]
    drop --> crown{"How old is the genus?<br/>(crown age)"}
    crown --> date["Stretch the tree so its<br/>root sits at that age and<br/>all tips at the present"]
```

A genus and its relatives are rarely separated perfectly by barcodes: a few BINs of the genus may fall among the related genera, or the other way round. The root is therefore put where the fewest BINs are on the wrong side; BINs of the genus that end up among the relatives stay in the tree as its earliest branches and are listed in the provenance, as some of them will be misidentified. When several places are equally good, the one that makes the genus tree most clock-like decides (minimal ancestor deviation, which is also used when there are no relatives to root with). Because the root is often uncertain, the provenance also says, per genus, how many other places are nearly as good, and how much the crown age and the genus's PD would change if the root were put there or if the genus's BINs among the relatives were left out.

The **crown age** is when the living members of the genus last shared an ancestor. It comes, in order of preference, from:

```mermaid
flowchart LR
    q{"Crown age<br/>known?"} -->|"from the user's<br/>literature table"| use1["use it"]
    q -->|"from TimeTree,<br/>and plausible"| use2["use it"]
    q -->|"otherwise"| rate["estimate it from<br/>how much the BINs differ,<br/>assuming a rate of change"]
    use1 --> cap["never older than<br/>90% of the stem age"]
    use2 --> cap
    rate --> cap
```

A TimeTree estimate counts as plausible only if it is younger than the stem age and not wildly older than what the barcodes suggest. Without that check, genera that are not natural groups (where TimeTree's species span a much older split) received crown ages of hundreds of millions of years.

## Step 4: grafting

```mermaid
flowchart LR
    subgraph before["Backbone"]
        p1(("split at<br/>S million<br/>years ago")) --- g["genus"]
        p1 --- s["sister genus"]
    end
    subgraph after["After grafting"]
        p2(("split at S")) ---|"stem: S minus C"| c(("crown at C"))
        c --- b1["BIN 1"]
        c --- b2["BIN 2"]
        c --- b3["BIN 3"]
        p2 --- s2["sister genus"]
    end
    before ==> after
```

The genus tip is replaced by a stem of length S − C (stem age minus crown age), followed by the genus tree, whose tips all end at the present. Every genus is handled independently, so the work runs in parallel and scales to tens of thousands of BINs.

## Step 5: samples

```mermaid
flowchart TB
    asv["A sequence from a sample"] --> search["Compare with all<br/>reference sequences"]
    search --> m{"How similar is<br/>the best match?"}
    m -->|"97% or more"| bin["It is that BIN"]
    m -->|"80 to 97%"| place["Add it to the tree<br/>next to its closest BIN,<br/>as far back in time as<br/>the difference suggests"]
    m -->|"below 80%"| out["Leave it out<br/>(reported)"]
    bin --> div["Diversity of the sample"]
    place --> div
```

With every sequence of a sample on the tree, diversity follows directly:

- **Faith's phylogenetic diversity**: the total length of the branches that connect a sample's organisms, in millions of years of evolutionary history.
- **Mean pairwise distance**: how long ago, on average, two organisms in the sample shared an ancestor.
- **UniFrac and PhyloSor**: how much of their branch length two samples share.

## Step 6: checking where barcodes stop being reliable

The method relies on barcodes being informative within genera. That can be checked: take pairs of BINs at increasing taxonomic distance (same genus, same family, same order), compare their barcode difference with their distance in a trusted tree built from other data, and see where the barcode difference stops growing. This is built into timegraft as a validation step, and it is the subject of planned work per animal phylum.

## What comes out

A dated tree with a BIN (or species hypothesis) at every tip, diversity values per sample and per pair of samples, and tables that record, for every genus and every branching point, where its age came from: TimeTree, a family-level stand-in, interpolation, the barcode rate, or a literature value. That record matters, because the ages are only as good as their sources; see [assumptions and limitations](limitations.md).
