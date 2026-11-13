# Configuration

All settings live in `config/config.yaml`, which runs the toy example as shipped. For your own runs, keep that file and put the settings that differ in a second file:

```bash
snakemake --cores 8 --configfile config/config.yaml my_run.yaml
```

Nested blocks are merged, so a second file only needs the keys it changes.

## Inputs

| Key | Meaning |
| --- | --- |
| `records` | Record table, tab-separated: `record_id`, `unit_id`, `genus`, `sequence`; optionally `family`, `order`, `class`, `phylum` (used by the validation and the genetic code) and `genetic_code` (an NCBI table number that overrides the lookup for that record) |
| `bcdm.path` | A BOLD data package (BCDM, TSV or Parquet). When set, it replaces `records` |
| `bcdm.marker`, `bcdm.taxa`, `bcdm.countries` | Filters applied to the BOLD package: marker code, `rank: [names]`, `country/ocean` values |
| `backbone` | A dated, ultrametric Newick tree in millions of years with one tip per genus. Ignored when `backbone_build.enabled` is true |
| `crown_ages` | Optional table of known crown ages (Ma), the place for literature values: `genus` and any of `crown_age` (a point age, used as given), `crown_age_min` and `crown_age_max` (bounds on the crown age); other columns, such as a reference, are ignored. See [crown ages and rates](#crown-ages-and-rates). A starting table is in `data/crown_ages/literature.tsv` |
| `samples` | Optional table `sample`, `unit_id`, `abundance`; enables diversity metrics |
| `asvs`, `asv_table` | Optional ASV sequences (FASTA) and table `sample`, `asv_id`, `reads`, instead of `samples` |
| `validation.reference_tree` | Optional trusted tree whose tips are unit ids; enables the validation step |

## Building the backbone

| Key | Default | Meaning |
| --- | --- | --- |
| `backbone_build.enabled` | `false` | Build the backbone instead of reading `backbone` |
| `backbone_build.source` | `static` | `static`: OToL release files; `api`: OToL web API |
| `backbone_build.synthesis`, `.taxonomy` | `opentree16.1`, `ott3.7.3` | Release versions |
| `backbone_build.context` | `Metazoa` | Taxon to match genera in, to avoid homonyms (`Fungi` for ITS) |
| `backbone_build.per_side` | 2 | Genera per child clade sent to TimeTree for each node |
| `backbone_build.species_per_genus` | 3 | Species ids sent along with each genus |
| `backbone_build.crown_ages_from_timetree` | `true` | Estimate genus crown ages from TimeTree species |
| `backbone_build.root_age` | | Root age in Ma, only needed if TimeTree cannot date the root |
| `backbone_build.interpolation` | `pure_birth` | Ages for nodes TimeTree does not date: `pure_birth`, by the number of genera in the full clade in the synthetic tree (static source only; otherwise `even`), or `even`, spaced evenly as in BLADJ ([results](../science/results/2026-11-interpolation.md)) |
| `backbone_build.synthesis_support` | `true` | Static source only: record per backbone node how many source trees of the synthesis support or conflict with it; downloads the synthesis `annotations.json` ([provenance](provenance.md)) |
| `backbone_build.upgrade_cache` | `true` | Ask TimeTree again for the evidence behind answers cached by versions before #47, keeping the cached ages ([provenance](provenance.md)) |
| `backbone_build.col.enabled` | `true` | Resolve genera in the Catalogue of Life and derive NCBI and OTT ids from it ([provenance](provenance.md#catalogue-of-life-as-identifier-hub), #48) |
| `backbone_build.col.release` | `316691` | ChecklistBank dataset key of the COL release (COL26.10 XR); an alias such as `3LXR` follows the latest release |
| `backbone_build.col.join` | `col` | `col`: send a genus to TimeTree with COL's NCBI id where it is one of the ids OTT holds for it, or OTT has none; `ott`: always OTT's |
| `backbone_build.col.species` | `false` | Also resolve every species name of the records (`backbone/col_species.tsv`) |
| `backbone_build.col.workers` | 4 | Concurrent requests to ChecklistBank |
| `backbone_build.ages.method` | `fit` | Resolving conflicting TimeTree ages: `fit` (all together, weighted by their evidence, ages decreasing towards the tips) or `precedence` (drop an age not younger than its dated ancestor; before #52) ([provenance](provenance.md#backbone-ages-from-a-weighted-fit)) |
| `backbone_build.ages.fit.sigma_floor` | 0.1 | Smallest uncertainty of a calibration (standard deviation of its log age) |
| `backbone_build.ages.fit.sigma_single` | 0.5 | Uncertainty of an age from a single study estimate |
| `backbone_build.ages.fit.higher_taxa_factor` | 2.0 | Uncertainty multiplier for ages through enclosing taxa |
| `backbone_build.ages.fit.taxonomy_only_factor` | 1.5 | Uncertainty multiplier for nodes no source tree of the synthesis supports |
| `backbone_build.ages.fit.huber_k` | 1.5 | Calibrations further than this many sigma from the fit are down-weighted |
| `backbone_build.ages.fit.outlier_z` | 3.0 | ... and reported as outliers beyond this |
| `backbone_build.ages.fit.min_gap` | 0.01 | Smallest step in age per edge between calibrated nodes, as a share of the older age |
| `backbone_build.anchoring.rule` | `majority` | Genera that are not a node of the synthesis: `majority` (smallest clade with more than `share` of their species) or `mrca` (MRCA of all; before #51) ([provenance](provenance.md#anchoring-genera-that-are-not-a-node-of-the-synthesis)) |
| `backbone_build.anchoring.share` | 0.5 | Share of the species the anchor clade must hold, more than; at least 0.5 and below 1 |
| `backbone_build.anchoring.min_species` | 3 | Fewer species in the synthesis: the MRCA |
| `backbone_build.anchoring.prefer_sampled` | `true` | Anchor on the species in the run's records rather than all species in OTT |
| `backbone_build.anchoring.min_sampled` | 2 | Sampled species needed in the synthesis to use them |
| `backbone_build.concordance.enabled` | `true` | Placement concordance: `backbone/concordance.tsv`; needed for re-anchoring ([provenance](provenance.md#placement-concordance), #49) |
| `backbone_build.concordance.flag_level` | `class` | A lineage or merged-id check fails when the lineages part at this rank or higher: `family`, `superfamily`, `order`, `class`, `phylum`, `kingdom` |
| `backbone_build.concordance.max_share` | 0.0 | The barcode-neighbour check fails when at most this share of the cross-genus hits is in a phylum of the backbone sisters |
| `backbone_build.concordance.min_hits` | 5 | Hits needed for the barcode-neighbour check |
| `backbone_build.concordance.min_flags` | 2 | Failed checks needed to flag a genus |
| `backbone_build.reanchor.mode` | `apply` | Move genera that placement concordance flags, before dating: `apply`; `report` (decide and record only); `off` ([provenance](provenance.md#re-anchoring-flagged-genera), #50) |
| `backbone_build.reanchor.min_severity` | `class` | Only move genera whose lineages part at this rank or higher |
| `backbone_build.reanchor.overrides` | | Curated table (`genus`, `action` `reanchor` or `keep`, `target`, `source`, `note`), e.g. `config/placement_overrides.tsv` |
| `backbone_build.cache` | `resources/api_cache.json` | Cache of TimeTree answers; keep it between runs |

## General settings

| Key | Default | Meaning |
| --- | --- | --- |
| `marker` | `coi` | Which block under `markers` to use |
| `outdir` | `results/toy` | Output directory |
| `n_shards` | 4 | Parallel jobs for the genus subtrees |
| `rooting.rule` | `tolerant` | `tolerant`: root on the edge with the fewest units on the wrong side; `strict`: the genus must be a clade with respect to the outgroup, else the fallback is used |
| `rooting.n_outgroup` | 10 | Outgroup units per genus. A top-level `n_outgroup` (the setting before #19) still overrides it |
| `rooting.outgroup_choice` | `spread` | `spread`: in turn from the genera of the nearest sister clades, longest sequences first; `sister`: the first units by id of the sister genus (before #19) |
| `rooting.outgroup_levels` | 3 | `spread`: how many sister clades to draw from, going up the backbone |
| `rooting.fallback` | `mad` | Without an outgroup, or when the strict rule fails: `mad` (minimal ancestor deviation) or `midpoint` |
| `rooting.uncertainty` | `true` | Record per genus how uncertain the root is and how much the crown age and PD depend on it (the `root_*`, `crown_age_min`/`_max` and `pd_rel_range` columns of `provenance.tsv`, [#34](https://github.com/naturalis/timegraft/issues/34)). Off leaves those columns empty and saves the extra dating of each alternative root |
| `max_records_per_unit` | 50 | Records considered when choosing a unit's representative |
| `genus_agreement` | 0.8 | Share of a unit's records that must agree on its genus |
| `prune_empty` | `true` | Drop backbone genera without units |
| `place_unassigned_units` | `true` | Place units without a usable genus next to their nearest unit |

Runs before #19 rooted with `rooting: {rule: strict, n_outgroup: 2, outgroup_choice: sister, fallback: midpoint}`. Why the defaults changed is in the [first results](../science/results/2026-11-iberia.md#why-rooting-fails).

## Marker plug-ins

Everything that differs between markers sits in `markers.<name>`:

| Key | COI default | ITS default | Meaning |
| --- | --- | --- | --- |
| `filter.min_length` | 200 | 100 | Minimum sequence length without gaps |
| `filter.max_ambiguous` | 0.01 | 0.01 | Maximum share of ambiguous bases |
| `filter.genetic_code` | 5 | | Default translation table for the reading-frame check (COI only); see [genetic codes](#genetic-codes) |
| `filter.genetic_codes` | `{}` | | Overrides, `rank: {name: code}` |
| `filter.genetic_code_table` | | | Your own name-to-code mapping instead of the shipped one |
| `align` | `none` | `none` | `none` (input already aligned) or `mafft` |
| `distance.method` | `k2p` | `kmer` | `k2p` needs aligned sequences; `kmer` takes `k` |
| `screening.enabled` | `true` | `false` | Set aside units that do not fit their genus before subtrees are built; see [screening units](#screening-units) |
| `screening.margin` | 0.02 | | A unit at least this much more similar (identity) to another genus than to any congener is set aside |
| `screening.min_identity` | 0.80 | | A unit with no congener at this identity is set aside, if the genus has units that do |
| `screening.query_cov` | 0.7 | | Share of a representative that must align in the search |
| `crown_age.method` | `rate` | `fraction` | How to set a genus crown age when no better source exists |
| `crown_age.rate` | 0.0115 | | Default substitution rate, per site per Myr, per lineage; sets the lower bound of the crown age, see [crown ages and rates](#crown-ages-and-rates) |
| `crown_age.rates` | `{}` | | Rate per taxon, `rank: {name: rate}` |
| `crown_age.rate_table` | | | The same as a table (`rank`, `name`, `rate`), as written by `timegraft.calibrate` |
| `crown_age.fraction` | | 0.5 | Crown age as a fraction of the stem age |
| `crown_age.cap_fraction` | 0.9 | 0.9 | A crown age never exceeds this fraction of the stem age |
| `crown_age.timetree_bound` | `max` | `max` | `max`: a TimeTree crown age is an upper bound on the crown age; `none`: TimeTree crown ages are not used |
| `crown_age.estimate` | `timetree_midpoint` | `timetree_midpoint` | The age between the bounds: `timetree_midpoint` (TimeTree where available, otherwise the midpoint, or for ITS the fraction of the stem age), `minimum` or `midpoint` (COI only); see [crown ages and rates](#crown-ages-and-rates) |
| `dating.method` | `mpl` | `mpl` | `mpl` (mean path length), `lsd2`, or `chronos` (penalised likelihood; see [dating within a genus](#dating-within-a-genus)) |
| `dating.seq_length` | 658 | 500 | Barcode length, for LSD2's variance model |
| `dating.model`, `dating.lambda` | | | `chronos` only: `correlated` (default), `relaxed`, `discrete` or `clock`, and the smoothing parameter (default 1) |
| `placement.match_identity` | 0.97 | 0.97 | At or above: an ASV counts as the unit it matches |
| `placement.place_identity` | 0.80 | 0.80 | At or above: attached near its nearest unit; below: unplaced |
| `placement.query_cov` | 0.9 | 0.8 | Share of the query that must align |
| `placement.height` | `rate` | `terminal` | How to set the attachment age |
| `placement.stem_fraction` | 0.5 | | Units set aside by screening that match nothing: joined to the stem of their own genus, this fraction of the way from crown to stem node |

### Screening units

Within-genus barcode distances set the lower bound of a genus's crown age, so a unit that does not belong in its genus (a misidentification, a contaminant, or a genus that is not a natural group) can make the genus look much older. Before genus subtrees are built, `screen_units` searches all representatives against each other with VSEARCH and sets a unit aside when:

- **`other_genus`**: it is at least `margin` more similar to a unit of another genus than to any unit of its own genus; or
- **`isolated`**: no unit of its own genus reaches `min_identity`, while other units of the genus do.

Genera with one unit are not screened, and a genus never loses all its units: if every unit would be set aside, all are kept (`kept_whole_genus`). Set-aside units are not discarded. They are placed like units without a usable genus: next to their nearest unit, at the age implied by their distance to it, or, if nothing reaches `placement.place_identity`, on the stem of their own genus (`stem_fraction` of the way from its crown to its stem node). They stay in the tree and in diversity, but do not shape the genus subtree or its crown age. `units/screened.tsv` gives every decision with its evidence.

The thresholds are judgements. On the Iberian data the search takes two minutes for 25,000 units, and these settings set aside about 5% of units ([results](../science/results/2026-11-screening.md)). Screening is off for ITS until the thresholds are tested on fungal data: identities across fungal genera are lower and less comparable.

### Crown ages and rates

The crown age of a genus lies between two bounds, and both are recorded ([#36](https://github.com/naturalis/timegraft/issues/36)):

1. A point age for the genus in the `crown_ages` table is used as it is (if it is not below the stem age, it is cut to `cap_fraction` of it).
2. Otherwise the **upper bound** is the smallest of `cap_fraction` of the stem age, the TimeTree crown age of the genus (`timetree_bound: max`), and `crown_age_max` of the table, which replaces TimeTree for that genus.
3. The **lower bound** (`method: rate`) is the age from the barcodes: the mean root-to-tip distance of the rooted genus tree divided by the substitution rate. Barcode distances saturate, also within genera, so with a rate from young calibrations this is a minimum, not an estimate. `crown_age_min` of the table raises it. With `method: fraction` there is no rate and no lower bound.
4. `estimate` picks the age between the bounds. `timetree_midpoint` (the default, [#42](https://github.com/naturalis/timegraft/issues/42)) takes the TimeTree crown age wherever it is the upper bound, that is, below `cap_fraction` of the stem age and with no table bounds for the genus, even when the barcode minimum is older; for the other genera it takes the midpoint, the geometric mean of the two bounds (or, with `method: fraction`, which has no lower bound, the fraction of the stem age). TimeTree aggregates published timetrees and is the best estimate there is where it covers a genus. The midpoint is the best estimate without it: it has the smallest worst-case relative error when only the bounds are known, and on the Iberian genera it was within a factor of two of TimeTree's crown age for three quarters of genera, the minimum for a quarter ([results](../science/results/2026-11-dating.md)). `minimum` takes the age from the barcodes, or `fraction` of the stem age, moved to the nearest bound it violates; `midpoint` uses the midpoint for every genus.

A TimeTree crown age is the MRCA of the species TimeTree knows under the genus name. It is far too old when the genus is not a natural group, which is why it only ever lowers an age. Earlier versions accepted it outright unless it exceeded `timetree_max_ratio` times the barcode age; that threshold made ages jump and is gone (the key is ignored, with a warning in the log).

`provenance.tsv` records the result. `crown_source` is `table`, `timetree`, `rate`, `fraction` or `midpoint`; `rate` and `fraction` carry a suffix when a bound was applied (`_capped` for the stem age, `_timetree_max`, `_table_max`, `_table_min`). `crown_lower` and `crown_upper` are the bounds, `crown_upper_source` says which source set the upper one (`stem`, `timetree` or `table`), and `barcode_age` is the age from the barcodes before any bound. The diversity step reports Faith's PD and MPD with every crown at its lower and at its upper bound next to the values for the tree itself, so the part of a result that rests on crown ages is visible ([outputs](running.md#outputs)).

The rate is resolved per genus from its taxonomy, from genus up to phylum (the most common family, order, class and phylum among the genus's records, written to `genera/taxonomy.tsv`): the first name found in `crown_age.rates` or, failing that, in `crown_age.rate_table` decides; otherwise `crown_age.rate`. `provenance.tsv` records the rate and where it came from (`rate`, `rate_source`). Use this for rates from the literature, for example a slower rate for elasmobranchs.

The default of 0.0115 is within the range of COI rates calibrated on splits of a few million years (0.007 to 0.026 per lineage), which is what a lower bound needs. `python -m timegraft.calibrate` derives rates per class from the genera of a run with a known crown age (a point age in the table or a TimeTree crown age below the stem age), as trimmed medians of root-to-tip distance over age. **Be careful with it**: on the Iberian genera the distance hardly grows with crown age, so the result is not a rate but the distance of a typical genus over the crown age of a typical genus, two to seven times slower than published rates, and it moved every genus towards the same age. Those rates are kept in `data/iberia/calibration/` as a record and are not used.

```bash
python -m timegraft.calibrate --provenance results/my_run/provenance.tsv \
    --taxonomy results/my_run/genera/taxonomy.tsv \
    --crown-ages data/crown_ages/literature.tsv --out my_rates.tsv
```

### Dating within a genus

Once the crown age is set, the nodes below it are dated by mean path length (`mpl`, built in) or LSD2 (`lsd2`, built from source), which both assume a clock within the genus. Penalised likelihood, which lets the rate vary over the tree, is available as `dating.method: chronos` (the `chronos` function of the R package ape; the workflow then uses `workflow/envs/chronos.yaml`, or R with ape on the path). It is not the default. On the Iberian genera it gave nearly the same within-genus branch lengths as mean path length, its own information criterion preferred the strict clock in every genus tested, and it cannot estimate the crown age: with all tips at the present and only bounds on the root, rate and time are not separable ([results](../science/results/2026-11-dating.md)). Genera with fewer than four units are dated by mean path length under any setting.

### Genetic codes

The COI reading-frame check needs the mitochondrial genetic code of each record, and a broad data set spans several (2 for vertebrates, 13 for tunicates, 9 for echinoderms and most flatworms, 4 for cnidarians and sponges, 5 for most other invertebrates). The code is resolved per record:

1. a `genetic_code` value in the record table, if there is one;
2. otherwise, going up from genus through subfamily, family, order, class, phylum and kingdom, the first name found in `filter.genetic_codes` or in the mapping (at the same rank, `filter.genetic_codes` wins);
3. otherwise `filter.genetic_code`.

The mapping ships as `src/timegraft/data/genetic_codes.tsv` (columns `name`, `rank`, `code`, `codes`, `ncbi_taxid`). Names are matched at any rank, because BOLD and NCBI do not always agree on rank. A name marked `mixed` covers more than one code (e.g. Chordata) and is skipped, so that a more specific name decides. The mapping holds every taxon from kingdom to order, plus families and genera wherever their code differs from what the group above would give (e.g. the families with table 14 inside Nematoda, which is mixed). To regenerate the mapping from the [NCBI taxonomy dump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz):

```bash
python -m timegraft.geneticcode taxdump.tar.gz src/timegraft/data/genetic_codes.tsv
```

`filtered.tsv` and `rejected.tsv` record the code used for each record and where it came from (`column`, `config`, `table` or `default`). The `filter_records` log lists the stop-codon rejection rate per code: a code with a much higher rate than the others suggests a wrong mapping.

## Validation

| Key | Default | Meaning |
| --- | --- | --- |
| `validation.ranks` | `[genus, family, order]` | Strata, from lowest to highest |
| `validation.per_stratum` | 200 | Pairs sampled per stratum |
| `validation.methods` | `k2p`, `kmer` | Distance methods to compare, each with its options |
