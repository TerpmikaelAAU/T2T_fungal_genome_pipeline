# T2T fungal genome assembly pipeline

> **This is a work in progress.**

## Article
Paper: "Gapless near Telomer-to-Telomer Assembly of Neurospora intermedia, Aspergillus oryzae, and Trichoderma asperellum from Nanopore Simplex Reads"
https://doi.org/10.3390/jof11100701

The tagged `v1.0.0` release on `main` is what the paper describes. This
branch (`v2-flexible-inputs`) generalises the pipeline beyond the paper's
~40 Mb fungal genomes: flexible entry points, a read-filtering grid, and
resource requests that scale with actual input size instead of fixed
per-rule numbers. Also tested end-to-end on a much larger, heavily
duplicated ~240 Mb genome (~75% transposable elements).

## What this pipeline does
Given raw Nanopore reads (or an already-basecalled BAM, or already-basecalled
FASTQ), the pipeline:
1. Optionally basecalls (pod5 entry point only) and trims adapters.
2. Filters reads across a grid of quality/length cutoffs, assembling each
   cell independently with `hifiasm --ont`.
3. Picks a winning assembly across the grid -- **twice, two independent
   ways**: fewest contigs, and highest BUSCO completeness (see "The
   read-filtering grid" below).
4. Polishes each winner with `dorado polish` where a BAM is available
   (once per distinct assembly, when several selectors pick the same one).
5. Runs BUSCO on each winner -- before and after polishing when it is
   polished -- and reports read/assembly stats and coverage.
6. Optionally recovers an organelle genome via `flye` + `GetOrganelle`.
7. Finds each sample's telomere motif (tidk), reports telomeres on every
   final assembly, optionally re-assembles with the motif
   (`hifiasm --telo-m`), and adds a third winner: most telomere-to-telomere
   contigs (see "Telomeres" below).
8. Screens every final assembly for contamination with Tiara (bacterial /
   archaeal contigs; see "Contamination screen" below).
9. Optionally annotates every final assembly the BAGS way: geneML genes,
   antiSMASH secondary-metabolite clusters, BUSCO on the proteins (see
   "Annotation (BAGS)" below).
10. Places every sample's final assemblies in a UFCG core-gene
   phylogeny alongside one NCBI genome per fungal genus (see "Phylogeny"
   below).

Everything is driven by Snakemake against the SLURM executor plugin -- one
job per rule instance, submitted and tracked by Snakemake itself.

## Requirements
- Nanopore reads, ideally kit 10.4 chemistry or newer.
- Reasonable coverage for your genome (the source paper used 100x+ at Q10;
  less can work, but sequence deep if you can).
- Access to the AAU BioCloud cluster (https://cmc-aau.github.io/biocloud-docs/),
  or another Snakemake 9 + SLURM-executor-plugin cluster with an equivalent
  profile.

## Installation

**0. Get the repo onto the cluster.**
Clone it, or download the zip via GitHub's "Code" button and upload/unzip it
on the server. Work from inside the repo directory for everything below.

**1. Create the Snakemake conda environment.**
```
conda env create -n snakemake_run -f Snakemake_env.yml
conda activate snakemake_run
```
This pins Snakemake 9.26.1 and the native SLURM executor plugin -- the
pipeline will not run correctly under Snakemake 7/8 or without the plugin.

**2. Dorado.**
Nothing to do by default: rule `get_dorado` fetches and unpacks a pinned
dorado release (`dorado.version` in `config.yaml`, default `2.1.2`) the first
time it's needed, and the result is `protected()` so it's only ever
downloaded once. To use an install you already have instead, set
`dorado.path` in `config.yaml`.

**3. Cluster profile.**
`profile/config.yaml` is already set up for BioCloud's SLURM executor
plugin. You shouldn't need to touch it unless adapting the pipeline to a
different cluster (see the comments in that file for what each setting
does).

## Try it first with example data
`config/config.yaml` ships pointing at `example_data/` -- tiny, fake reads
that let you run the whole pipeline (containers, scheduling, all the rules)
without your own data first. See `example_data/README.md`. Swap in your own
sample paths when you're ready for a real run.

## Configuring a run

All configuration lives in `config/config.yaml`; it's commented in place, but
the shape is:

```yaml
samples:
  my_sample:
    type: fastq            # pod5 | bam | fastq -- see entry points below
    path: /abs/path/to/data   # or a path relative to the repo root
    busco_lineage: fungi_odb12   # optional, default shown
    organelle: false              # optional, recover an organelle genome
    organelle_type: fungus_mt     # REQUIRED if organelle: true -- GetOrganelle's -F type
    dorado_correct: false         # optional, per-sample override of dorado_correct.enabled below
    ultralong: false               # optional, per-sample override of ultralong.enabled below
    subsample: {coverage: 100, genome_size: 40mb}  # optional coverage cap
    porechop: true                 # optional, adapter trimming (default on)
    filter:                        # optional, PER-SAMPLE grid override
      min_q:   [20]
      min_len: [1000]

dorado:
  version: "2.1.2"          # or: path: /my/own/dorado/bin/dorado

dorado_correct:
  enabled: false             # hifiasm --ont does its own correction; see below
  min_q: 10                  # fixed cutoff fed INTO dorado correct
  min_len: 10000              # (the filter grid's min_len still applies AFTER)
  index_size: "4G"
  inference_device: "cuda:all"

ultralong:
  enabled: false             # hifiasm --ul, a second fixed-cutoff read set
  min_q: 10
  min_len: 50000

hifiasm_min_input_mb: 1      # below this, skip hifiasm instead of SIGILL-ing

filter:                      # the DEFAULT (min_q, min_len) grid
  min_q:   [10, 15, 20]
  min_len: [1000, 2500, 7500, 10000]
```

### Entry points (`type:`)
| `type:` | `path:` must be | Runs |
|---|---|---|
| `pod5`  | a directory of `.pod5` files | basecall -> filter -> assemble -> polish |
| `bam`   | one basecalled `.bam` file | filter -> assemble -> polish |
| `fastq` | one `.fastq` or `.fastq.gz` file | filter -> assemble (no polishing -- there's no move table to polish with) |

### Duplex data
Duplex reads from `dorado duplex` work with the `bam` and `fastq` entry
points, with no extra config. Such a file holds every duplex read *and* the
two simplex reads it was made from, so the pipeline handles it like this:
- **Assembly:** the simplex parents of duplex reads are dropped, so each
  molecule enters hifiasm once -- by the `dx:i:-1` tag for a BAM, and by
  name for a FASTQ (duplex reads are named `<id1>;<id2>`, so `id1` and
  `id2` go). Simplex reads without a duplex partner are kept. The numbers
  are in `logs/bam_to_fastq/<sample>.log` or `logs/prepare_fastq/<sample>.log`.
  For a FASTQ the read names must be dorado's own; renamed reads can't be
  matched to their parents.
- **Polishing:** `dorado polish` rejects duplex data, so it gets the simplex
  reads only (parents included -- they cover the same molecules) and the
  duplex read group is removed from the aligned BAM's header.

Simplex data passes through both steps unchanged.

### The read-filtering grid, and the per-sample override
`chopper` filters reads at every `(min_q, min_len)` combination in `filter:`,
and each cell gets assembled independently with `hifiasm`. The top-level
`filter:` is the default for every sample that doesn't set its own (see
`sample_a` in `config/config.yaml`, which overrides it to a single cell).
**Give a sample its own `filter:` block if it needs a different grid than
the others** -- otherwise editing the shared default changes what every
*other* sample's `contig_count`/`highest_busco`/`assembly_stats` requires
too, including grid cells that sample never ran.

For a sample with `dorado_correct` on, this grid's `min_q` is ignored --
only its `min_len` values apply, and only AFTER correction -- see below.

**Two independent rules pick a winner across that grid, in parallel:**
- `contig_count` picks the fewest-contigs assembly. Check the winner's total
  length in `assembly_stats.tsv` against your expected genome size -- for a
  repeat-rich genome, fewer contigs can mean collapsed repeats rather than a
  genuinely better assembly.
- `highest_busco` runs BUSCO on *every* grid cell's assembly and picks the
  highest-completeness one instead -- doesn't share `contig_count`'s
  collapsed-repeat blind spot. **Cost warning:** this is one BUSCO run per
  grid cell (e.g. 12 for a 3x4 grid), on top of the final BUSCO run each
  winner gets afterwards (two for pod5/bam samples: before and after
  polishing).

Both winners go through the rest of the pipeline independently -- polishing
if a BAM is available (with BUSCO run both before and after, so
`summary.txt` shows what polishing changed), a final BUSCO run, stats, and a
deliverable -- so
every sample ends up with two directly comparable results:
`results/<sample>/lowest_contig/` and `results/<sample>/highest_busco/`.
Neither is automatically "more correct"; compare their `summary.txt` files
(contig count, N50, total length vs. expected genome size, and BUSCO score)
and use your judgement.

Grid cells smaller than `min_assembly_mb` (config, default 10 Mb) are never
picked by any selector, and BUSCO is skipped on them (scored 0%): on sparse
data a single short contig would otherwise "win" fewest contigs, and BUSCO
can crash on such a gene-less assembly. If no cell reaches it, the selectors
fall back to every non-empty cell and say so in `logs/select/`. Each
selector's log there lists every cell's contigs and size.

### `dorado_correct` vs hifiasm's own correction
`hifiasm --ont` does its own ONT-specific read correction, and that's the
default (`dorado_correct.enabled: false`, overridable per sample with the
`dorado_correct:` key). Only turn it on to experiment with dorado's
correction instead: on a low-N50 (~5 kb) test dataset it discarded ~85% of
reads, leaving too little coverage for hifiasm to assemble anything useful.
It may do better on data with a longer N50.

**Does the `min_q`/`min_len` grid still apply when this is on?** Half of
it. dorado correct itself always runs ONCE per sample, on the fixed
`dorado_correct.min_q`/`min_len` cutoff above (same idea as `ultralong`'s
fixed cutoff) -- not the `filter:` grid. Its output is FASTA with no
per-base quality (correction discards it), so it's padded with a
placeholder quality string just so hifiasm can read it as FASTQ. *That*
padded FASTQ is then filtered again by the `filter:` grid's `min_len`
values (its own `filter:` override if the sample has one) -- one assembly
per length, same as any other sample's grid -- but not by `min_q`, since
there's no real quality left to filter on.

### BUSCO lineage
Set `busco_lineage` per sample to match the organism -- e.g. a
`stramenopiles_odb*` dataset for an oomycete like *Phytophthora*, which is
not a fungus despite this pipeline's name. `busco --list-datasets` (inside
the BUSCO container -- see `workflow/rules/7_BUSCO.smk` for the image) lists
what's available. Don't compare BUSCO scores across dataset versions
(odb10 vs odb12).

### Telomeres
With `telomere: enabled: true` (the default) -- see
`workflow/rules/5_1_telomere_motif.smk`, `6_02_most_t2t.smk` and
`9_1_telomere_report.smk`:

1. **Find the motif, for any genome.** `tidk explore` runs on the sample's
   first-pass assembly with the fewest contigs and lists candidate repeat
   units near contig ends. The top `top_motifs` (10) are each checked with
   `tidk search`; the motif that makes the most contig ends telomeric wins.
   The ranking is in `results/<sample>/telomere_motif_candidates.tsv`.
2. **Re-assemble with it (optional).** If that motif covers at least
   `rerun_min_fraction` (80%) of all telomeric ends any candidate found, and
   at least `rerun_min_ends` ends, every grid cell is re-assembled with
   `hifiasm --telo-m <motif>`. The selectors then choose across BOTH
   hifiasm runs (`+telo` rows in `summary.txt`). This doubles the hifiasm
   cost when it fires; set `rerun_hifiasm: false` to skip it.
3. **Third selector, `most_t2t`.** The grid cell with the most
   telomere-to-telomere contigs (ties: more telomeric ends, then fewer
   contigs) -- polished, BUSCO'd and reported like the other two, under
   `results/<sample>/most_t2t/`.
4. **Report on every final assembly:** `telomeres.tsv` (per contig:
   repeats at each end, telomeric?, T2T?), a `tidk plot` SVG under
   `telomeres/`, and a TELOMERES section in `summary.txt`.

An end counts as telomeric with at least `min_repeats` (10) motif copies in
its terminal `window` (2000 bp). If no candidate finds any telomere,
`fallback_motif` (TTAGGG) is used for the counts and hifiasm is not re-run.

### Phylogeny
With `phylogeny: enabled: true` (the default), the pipeline also builds one
tree over all samples: both final assemblies per sample (`lowest_contig` and
`highest_busco`) plus **one NCBI genome for every genus** under
`phylogeny: taxon:` (default `Fungi`). It's a protein tree from
[UFCG](https://ufcg.steineggerlab.com)'s universal fungal core genes (the
tool used for the paper, see `scripts/UFCG_Phylogeny/`): `ufcg profile`
per genome, `ufcg align` (MAFFT) + concatenation, then FastTree (LG+gamma).

- The NCBI download (`workflow/rules/0_4_ncbi_genus_representatives.smk`)
  and the UFCG profiling of those reference genomes
  (`workflow/rules/10_phylogeny.smk`) don't depend on any sample, so they
  start right away and run in parallel with basecalling/assembly. Only the
  samples' own profiles and the tree wait for the assemblies.
- Genus comes from NCBI taxonomy. Per genus, the representative is the NCBI
  "reference genome" if there is one, then the best assembly level
  (Complete > Chromosome > Scaffold > Contig), then the highest contig N50.
  The picks are listed in `data/phylogeny/genus_representatives.tsv`.
- **Cost:** all of Fungi is a few thousand genomes, **~100 GB** to
  download, and one UFCG profile job per genome. Everything under
  `data/phylogeny/` is kept so this only ever happens once; delete that
  folder to start over (e.g. after changing `taxon`).
- The reference profile jobs never hold up the assemblies: they rank below
  every other job, and at most `ufcg_reference_slots` (10, in
  `profile/config.yaml`) of the profile's 30 SLURM job slots run them at
  once.
- A genome UFCG can't profile is left out of the tree with a warning (see
  `logs/phylogeny/`) rather than failing the run.
- **Skipped automatically** when every sample is bundled `example_data/`.
- Optional: `export NCBI_API_KEY=<key>` before running Snakemake; the
  `datasets` CLI picks it up and NCBI then allows more requests per second.
- One genome per genus places a sample at genus/family level; it is not
  meant to resolve species. UFCG's markers are fungal, so this doesn't
  suit oomycetes.

### Contamination screen
With `contamination: enabled: true` (the default), every final assembly is
run through [Tiara](https://github.com/ibe-uw/tiara)
(`workflow/rules/8_contamination_tiara.smk`), a small deep-learning
classifier that labels each contig as eukarya, bacteria, archaea,
organelle (mitochondrion/plastid) or unknown. `summary.txt` gets a
CONTAMINATION section with contigs and bp per class, and a warning listing
any bacterial/archaeal contigs; per-contig classes and probabilities are in
`results/<sample>/<selector>/contamination/tiara.tsv`. It only reports --
nothing is removed.

Why Tiara and not NCBI's FCS-GX: FCS-GX needs a ~500 GB database held in
RAM to be fast, and NCBI runs it on every genome submission anyway. Tiara
ships its models, needs no database and takes minutes. Contigs shorter than
`min_len` (3 kb) aren't classified. Tiara isn't on bioconda and pins old
dependencies (Python ≤ 3.9), so the first run installs it once into
`resources/tiara-<version>/` inside the official `python:3.9-slim` image.

### Annotation (BAGS)
Off by default; set `bags: enabled: true`. The steps of the
[BAGS pipeline](https://github.com/TerpmikaelAAU/BAGS), run on every
selector's final assembly (`workflow/rules/11_annotation_bags.smk`):
geneML gene prediction → antiSMASH with geneML's genes
(`--genefinding-tool none`, plus `--cc-mibig --cb-general`) → BUSCO in
protein mode, which measures how complete the *annotation* is (the
genome-mode BUSCO measures the assembly). Output per assembly in
`results/<sample>/<selector>/annotation/` (GFF3, proteins, antiSMASH
folder, protein BUSCO), and one row per assembly in
`results/annotation_summary.tsv`: genes, BUSCO, and BGC regions by type
(NRPS, PKS, terpene, RiPP, other; a hybrid region counts in each of its
classes).

geneML is only distributed on PyPI, so the first run installs it once into
`resources/geneml-<version>/` (like dorado) inside the official
`python:3.12-slim` image; the antiSMASH databases are likewise downloaded
once into `resources/antismash_databases/`.

## Containers
Every tool (BUSCO, hifiasm, Flye, ...) runs inside a container, downloaded
automatically the first time that tool is needed and reused after that. You
don't need to install these tools yourself, and there is no setup step.

Each container is one pinned image from
[biocontainers](https://biocontainers.pro/) (e.g.
`quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0`), listed under
`container:` in each file under `workflow/rules/`. The one exception is
NCBI `datasets`, which is no longer on bioconda (StaPH-B's
`staphb/ncbi-datasets` image instead), and UFCG, whose official
`endix1029/ufcg` image ships the core gene database the bioconda one lacks,
geneML, installed from PyPI into `python:3.12-slim` (see "Annotation
(BAGS)"), and Tiara, likewise into `python:3.9-slim` (see "Contamination
screen"). `workflow/envs/*.yml`
still documents the same tool+version as a conda environment, but is no
longer used to install anything -- it's just a reference.

## Running the pipeline

**Run the workflow.** `slurm_submit.sbatch` wraps this in a single
lightweight SLURM job (1 CPU, 4G -- each *rule* still gets its own
appropriately-sized job via the executor plugin):
```
sbatch slurm_submit.sbatch
squeue --me      # watch it
```
Alternatively, since the orchestrator itself is cheap and just submits/polls,
run it directly on the login node inside `tmux` -- the more common way to
drive a long multi-day run interactively:
```
tmux new -s my_run
conda activate snakemake_run
snakemake --profile profile
# detach with Ctrl-b d; reattach later with: tmux attach -t my_run
```

**Multiple samples configured at once:** `rule all` builds every sample in
`config.yaml` in full. While iterating on one sample, target its outputs
explicitly instead, e.g.:
```
snakemake --profile profile \
  results/my_sample/lowest_contig/summary.txt \
  results/my_sample/highest_busco/summary.txt
```

**If a run dies and leaves the working directory locked:**
```
snakemake --profile profile --unlock
```

**Debugging a failed rule:** most intermediates are `temp()` and get deleted
as soon as their consumers finish, taking the evidence with them. Re-run
with `--notemp` to keep everything while debugging.

## Output

```
results/<sample>/
  read_stats.tsv                      # seqkit stats for the raw (and trimmed) reads -- shared, sample-level
  polish_groups.tsv                   # pod5/bam only: which selectors picked the same assembly (polished once)

  lowest_contig/                      # winner selected by fewest contigs
    summary.txt                       #   human-readable: reads, grid, coverage, BUSCO (before + after polishing if polished)
    assembly_stats.tsv                #   seqkit stats for every grid candidate + this winner
    <sample>_lowest_contig_final.fasta  # THE deliverable for this selector -- polished if a BAM was available

  highest_busco/                      # winner selected by highest BUSCO completeness
    summary.txt
    assembly_stats.tsv
    <sample>_highest_busco_final.fasta  # THE deliverable for this selector

  most_t2t/                           # only with telomere: enabled -- winner by most T2T contigs
    ...                               #   same files as above
  <selector>/telomeres.tsv            # per-contig telomeres (telomere: enabled), + telomeres/*.svg
  telomere_motif_candidates.tsv       # tidk explore candidates and the chosen motif

  <selector>/contamination/tiara.tsv  # Tiara class per contig (contamination: enabled)
  <selector>/annotation/              # only with bags: enabled -- GFF3, proteins, antismash/, protein BUSCO

results/annotation_summary.tsv        # only with bags: enabled -- genes, BUSCO, BGCs per final assembly

results/phylogeny/                    # only if phylogeny: enabled: true
  genus_tree.nwk                      # Newick tree: all samples' final assemblies + one NCBI genome per genus
  ufcg_concatenated_alignment.fasta   # the concatenated UFCG core-gene protein alignment behind it
  tip_labels.tsv                      # tip label -> accession, organism, phylum/class/order/family (e.g. for iTOL)

logs/                   # per-rule, per-sample logs
snake_log/               # the Snakemake orchestrator's own run logs
```

All selectors run all the way through independently -- there is no single
"the" final assembly, there are two (three with telomeres on): compare
`lowest_contig/`, `highest_busco/` and `most_t2t/summary.txt` and pick
whichever looks better for your genome.

`data/` holds intermediates and is mostly cleaned up automatically once
nothing downstream needs it. These persist deliberately:
`data/dorado_basecall/<sample>.bam` (basecalling is expensive to redo) and,
if a sample sets `organelle: true`, `data/getorganelle/<sample>/Mitochondria`
(the only copy of that result). With the phylogeny on,
`data/phylogeny/` is kept too, so the ~100 GB NCBI download and the
reference genomes' UFCG profiles are only made once.

## Repo layout
- `Snakefile` -- entry-point handling, input resolvers, adaptive resource
  helpers, and the rule includes/targets. Commented throughout; read it
  before adding a new rule or entry point.
- `workflow/rules/*.smk` -- one file per pipeline stage, numbered in
  roughly the order they run in the DAG. Each has a header comment
  explaining what it does and when it's reached.
- `workflow/envs/*.yml` -- one conda env per tool. Kept as reference only;
  see "Containers" above for what actually runs each tool.
- `workflow/scripts/summarize.py` -- builds
  `results/<sample>/<selector>/summary.txt`.
- `workflow/scripts/drop_duplex_parents.sh` -- copies a FASTQ sample's
  reads without the simplex parents of duplex reads (see "Duplex data").
- `workflow/scripts/pick_highest_busco.py` -- picks the `highest_busco`
  selector's grid winner (see `rule highest_busco`).
- `workflow/scripts/telomeres.py` (+ `pick_telomere_motif.py`,
  `pick_most_t2t.py`, `telomere_report.py`) -- reads `tidk search` output:
  motif choice, the `most_t2t` pick, and the per-contig telomere table.
- `workflow/scripts/bags_summary.py` -- one annotation row per assembly
  (see "Annotation (BAGS)").
- `workflow/scripts/pick_genus_representatives.py` -- picks one NCBI
  assembly per genus for the phylogeny (see `rule pick_genus_representatives`).
- `config/config.yaml`, `profile/config.yaml` -- see above.
- `example_data/` -- tiny fake reads for a first test run; see its README.
- `scripts/` -- separate, standalone scripts used to make the paper's
  figures (phylogeny, circos plots, telomere/mitochondria checks, ...).
  Barebones and not part of the Snakemake DAG; see `scripts/README.md`.

## Automatic checks
`.github/workflows/dry-run.yml` runs on every push and pull request: it
compiles the Python scripts and runs `snakemake -n` on the default config,
with every optional feature on, and with every feature off. It runs no
tools and needs no cluster -- it catches a rule or config change that breaks
the workflow before it's merged.

## Cluster notes (AAU BioCloud)
Learned from running this pipeline on BioCloud; useful if you're adapting
`profile/config.yaml`/`slurm_submit.sbatch` to a different cluster too.

- **`--partition` has no effect here.** BioCloud assigns partitions
  automatically from each job's memory-per-CPU ratio and node features.
- **GPU:** exactly one node, `bio-node10` -- 1x NVIDIA A10 (Ampere, CC 8.6,
  23 GB VRAM), 64 threads, 256 GB RAM. Request it with `--gres=gpu:a10:1`
  (BioCloud's GPU_GRES in the Snakefile), not via `--partition`.
- **Interactive jobs land on the `interactive` partition** (`bio-node11`, no
  GPU) unless you explicitly pass `--gres` for GPU work interactively.
- **Fat CPU nodes** (`zen5x` ~2.3 TB RAM/288 threads, `zen3x` ~2 TB) are
  better than the GPU node for memory-hungry non-GPU work like hifiasm or
  dorado correct's overlap stage.
- **Compute nodes have internet access** -- container downloads and dorado's
  own download both work from inside a job, not just the login node.
- **The SLURM account auto-guess can fail** (`sacct: invalid option -- '1'`)
  without breaking anything. If jobs start getting rejected for account
  reasons, set an explicit `slurm_account` under `default-resources` in
  `profile/config.yaml`.
- **Run the Snakemake orchestrator on the login node, not inside a SLURM
  job.** Snakemake warns about this; the orchestrator barely uses any
  resources itself (it only submits jobs and polls), so wrapping it in its
  own allocation has no benefit and can make that job look hung.

## Known gotchas
- **`temp()` deletes your debugging evidence.** Re-run with
  `snakemake --profile profile --notemp` while debugging a failing rule so
  its intermediates stick around for inspection.
- **A missing declared input crashes DAG construction**, even for a target
  that's already complete and up to date. If you delete a sample's raw data
  once its `results/` are final, a later run touching that sample can crash
  with `MissingInputException` trying to plan how to regenerate it. Fix:
  leave a 0-byte placeholder at the deleted path with an old timestamp
  (`touch -d "<old date>" <deleted path>`).
- **`--until <rule>`** stops the pipeline at a chosen rule, useful for
  building just far enough to debug one stage.
- **Changing a rule's `container:` tag downloads a new image** the next time
  that rule runs. Old ones accumulate under `.snakemake/singularity/`;
  `rm -rf .snakemake/singularity` after a round of version bumps to free the
  space.
- **`script:` directive paths resolve at DAG-build time** -- a typo'd or
  missing script path kills even a dry run (`-n`), before any rule executes.
- **Wildcard constraints matter for parsing.** `minq`, `minlen`, and `block`
  are constrained to `\d+` in the Snakefile so a filename like
  `sample_a_q10_l5000` parses as `input=sample_a, minq=10, minlen=5000`
  instead of Snakemake trying every possible split.
- **Every rule's `log:` directive must contain every one of that rule's
  wildcards**, or Snakemake refuses to build the DAG at all.
- **`git add -A` sweeps in generated junk.** `data/`, `results/`, `logs/`,
  etc. are gitignored, but always check `git status` before a broad `add`
  in case something generated ended up somewhere unexpected.
- Do NOT set your conda channel priority to `strict` -- this only matters if
  you fall back to `workflow/envs/*.yml` instead of containers, but strict
  priority makes those environments fail to resolve.

## How to cite
If this pipeline (or its output) contributes to a publication, please cite
the paper it comes from:

> Terp, M.; Nyitrai, M.; Rusbjerg-Weberskov, C.E.; Sondergaard, T.E.;
> Lübeck, M. Gapless near Telomer-to-Telomer Assembly of *Neurospora
> intermedia*, *Aspergillus oryzae*, and *Trichoderma asperellum* from
> Nanopore Simplex Reads. *Journal of Fungi* **2025**, *11*, 701.
> https://doi.org/10.3390/jof11100701

```bibtex
@article{terp2025gapless,
  author  = {Terp, Mikael and Nyitrai, Mark and Rusbjerg-Weberskov, Christian Enrico and Sondergaard, Teis E. and L{\"u}beck, Mette},
  title   = {Gapless near Telomer-to-Telomer Assembly of Neurospora intermedia, Aspergillus oryzae, and Trichoderma asperellum from Nanopore Simplex Reads},
  journal = {Journal of Fungi},
  year    = {2025},
  volume  = {11},
  number  = {10},
  pages   = {701},
  doi     = {10.3390/jof11100701}
}
```

If you're specifically using the `v2-flexible-inputs` generalisation (flexible
entry points, the read-filtering grid, adaptive resources) rather than the
exact paper pipeline, please also note the pipeline version/tag you ran
(e.g. `v2.0.0`) and link to this repository.
