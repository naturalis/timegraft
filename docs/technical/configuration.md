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
| `crown_ages` | Optional table `genus`, `crown_age` (Ma). Always used as given; the place for literature values |
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
| `crown_age.method` | `rate` | `fraction` | How to set a genus crown age when no better source exists |
| `crown_age.rate` | 0.0115 | | Substitutions per site per Myr, per lineage |
| `crown_age.fraction` | | 0.5 | Crown age as a fraction of the stem age |
| `crown_age.cap_fraction` | 0.9 | 0.9 | A crown age never exceeds this fraction of the stem age |
| `crown_age.timetree_max_ratio` | 3 | | Reject a TimeTree crown age above this many times the rate-based one |
| `dating.method` | `mpl` | `mpl` | `mpl` (mean path length) or `lsd2` |
| `dating.seq_length` | 658 | 500 | Barcode length, for LSD2's variance model |
| `placement.match_identity` | 0.97 | 0.97 | At or above: an ASV counts as the unit it matches |
| `placement.place_identity` | 0.80 | 0.80 | At or above: attached near its nearest unit; below: unplaced |
| `placement.query_cov` | 0.9 | 0.8 | Share of the query that must align |
| `placement.height` | `rate` | `terminal` | How to set the attachment age |

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
