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

Some BINs do not fit the genus they are assigned to: the records were misidentified, the sample was contaminated, or the genus is not a natural group. Because the differences between the BINs of a genus set how old the genus is taken to be (step 3), one such BIN can make a genus look much older than it is. Every BIN is therefore compared with every other, and a BIN is set aside when it is clearly more similar to a BIN of another genus than to any of its own genus, or when it is unlike all the others of its genus. A set-aside BIN is not thrown away: it is attached to the tree next to its closest relative, or, if it resembles nothing, on the branch leading to its own genus, so that it still counts in diversity but no longer sets the age of the genus.

## Step 2: the dated backbone

The backbone is a tree of genera with known ages:

1. All genera in the data are looked up in the Open Tree of Life, which gives how they are related.
2. For each branching point, TimeTree is asked when the lineages on either side split. TimeTree does not know every genus, so it is also asked about species within the genera and, as a last resort, about the families or orders they belong to.
3. Branching points that TimeTree cannot date are placed between the nearest dated points above and below by a pure-birth model: the more genera a clade holds in the Open Tree of Life, the older its expected age. The count includes genera that are not in the data, so the age of a clade does not depend on which other genera happen to be sampled ([results](results/2026-11-interpolation.md)). When those counts are not available, the points are spread evenly instead.

Each genus in the backbone then has a **stem age**: when it split from its closest relatives.

## Step 3: a small tree per genus

```mermaid
flowchart TB
    bins["The genus's BINs<br/>+ ten BINs from<br/>the closest related genera"] --> distances["How different is each pair?<br/>(barcode distances)"]
    distances --> join["Join the most similar<br/>step by step into a tree"]
    join --> root["Find the root:<br/>where the related genera<br/>attach"]
    root --> drop["Remove the related genera"]
    drop --> crown["How old is the genus?<br/>(crown age, see below)"]
    crown --> date["Stretch the tree so its<br/>root sits at that age and<br/>all tips at the present"]
```

A genus and its relatives are rarely separated perfectly by barcodes: a few BINs of the genus may fall among the related genera, or the other way round. The root is therefore put where the fewest BINs are on the wrong side; BINs of the genus that end up among the relatives stay in the tree as its earliest branches and are listed in the provenance, as some of them will be misidentified. When several places are equally good, the one that makes the genus tree most clock-like decides (minimal ancestor deviation, which is also used when there are no relatives to root with). Because the root is often uncertain, the provenance also says, per genus, how many other places are nearly as good, and how much the crown age and the genus's PD would change if the root were put there or if the genus's BINs among the relatives were left out.

The **crown age** is when the living members of the genus last shared an ancestor. Barcodes alone cannot give it. They say how much the BINs of a genus differ, and that difference stops growing after a while (saturation), also within genera: in half of the Iberian genera, the two most different BINs of the genus differ as much as the genus differs from its relatives. What the sources can give is a range:

```mermaid
flowchart LR
    q{"Crown age in the<br/>literature table?"} -->|"yes"| use1["use it"]
    q -->|"no"| low["at least: how much the BINs<br/>differ, at the rate of change<br/>measured on young splits"]
    q -->|"no"| up["at most: 90% of the stem age,<br/>and TimeTree's age for the genus"]
    low --> pick["choose within the range:<br/>TimeTree's age if it has one,<br/>otherwise the midpoint"]
    up --> pick
```

- **The barcodes give a minimum.** Dividing the difference between BINs by the rate at which COI changes over a few million years gives the age the genus must at least have. Saturation hides older change, so the true age can be much higher but not lower.
- **TimeTree's age for a genus is a maximum, and where it is below the stem, the answer.** TimeTree combines the dated trees published for a group, which is the best information there is. Its age for a genus is the age of the common ancestor of all species it knows under the genus name. For a natural group that is about right; for a genus that is not a natural group it spans a much older split, and can be hundreds of millions of years. It is therefore never allowed to raise a crown age above the stem age, but below that it is used as the crown age.
- **Where TimeTree has no age, the crown age is the middle of the range.** For about three quarters of genera TimeTree has no usable age. Their crown age is then the midpoint of the range on a logarithmic scale, between the barcode minimum and 90% of the stem age. Tested against the genera TimeTree does date, this midpoint was within a factor of two for three quarters of them; the barcode minimum for only a quarter.
- **The range is kept.** Both limits are recorded for every genus, and diversity is also computed with every genus at its minimum and at its maximum, so that the effect of this uncertainty on a result can be read off.
- **No thresholds.** An earlier version used TimeTree's age outright unless it was more than three times the barcode estimate. A small change in the data could then move a crown age by a factor of three. With limits, a small change in the data gives a small change in the age.

Below the crown, the tree is stretched so that its root sits at the crown age and all tips at the present, assuming that the BINs of a genus change at about the same rate (mean path length, or LSD2). Methods that let the rate vary over the tree (penalised likelihood) were tested and are available, but are not the default: within a genus they gave nearly the same branch lengths, and they cannot estimate the crown age itself. The [dating results](results/2026-11-dating.md) give the evidence.

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

A dated tree with a BIN (or species hypothesis) at every tip, diversity values per sample and per pair of samples, and tables that record, for every genus and every branching point, where its age came from: TimeTree, a family-level stand-in, interpolation, the barcodes with the rate that was used and any limit that was applied, or a literature value. That record matters, because the ages are only as good as their sources; see [assumptions and limitations](limitations.md).
