# T2T fungal genome assembly pipeline

## Article
Paper: "Gapless near Telomer-to-Telomer Assembly of Neurospora intermedia, Aspergillus oryzae, and Trichoderma asperellum from Nanopore Simplex Reads"
https://doi.org/10.3390/jof11100701

The tagged `v1.0.0` release on `main` is what the paper describes. This
branch (`v2-flexible-inputs`) generalises the pipeline beyond the paper's
~40 Mb fungal genomes -- flexible entry points, a read-filtering grid, and
resource requests that scale off actual input size (`scaled_mem`/
`scaled_time` in the Snakefile) instead of fixed per-rule numbers. That
scaling has been exercised end-to-end on a much larger, heavily duplicated
genome (~240 Mb, ~75% transposable elements -- typical of many oomycete
genomes), confirming memory/runtime requests hold up well outside the
paper's ~40 Mb fungal genomes, where the per-rule figures the paper's
pipeline shipped with were tuned.

## What this pipeline does
Given raw Nanopore reads (or an already-basecalled BAM, or already-basecalled
FASTQ), the pipeline:
1. Optionally basecalls (pod5 entry point only) and trims adapters.
2. Filters reads across a grid of quality/length cutoffs, assembling each
   cell independently with `hifiasm --ont`.
3. Picks a winning assembly across the grid -- **twice, two independent
   ways**: fewest contigs, and highest BUSCO completeness (see "The
   read-filtering grid" below).
4. Polishes each winner with `dorado polish` where a BAM is available.
5. Runs BUSCO on each winner and reports read/assembly stats and coverage.
6. Optionally recovers an organelle genome via `flye` + `GetOrganelle`.

Everything is driven by Snakemake against the SLURM executor plugin -- one
job per rule instance, submitted and tracked by Snakemake itself.

## Requirements
- Nanopore reads, ideally kit 10.4 chemistry or newer.
- Reasonable coverage for your genome (the source paper used 100x+ at Q10;
  less can work, as shown there, but sequence deep if you can).
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
downloaded once. If you'd rather use an install you already have, set
`dorado.path` in `config.yaml` (see below) and this rule never runs.

**3. Cluster profile.**
`profile/config.yaml` is already set up for BioCloud's SLURM executor
plugin -- unlike the paper-era pipeline, there is no partition to edit by
hand in every rule; BioCloud assigns partitions automatically from each
job's memory-per-CPU ratio. You generally shouldn't need to touch this file
unless you're adapting the pipeline to a different cluster (see the
comments in that file for what each setting does).

## Configuring a run

All configuration lives in `config/config.yaml`; it's commented in place, but
the shape is:

```yaml
samples:
  my_sample:
    type: fastq            # pod5 | bam | fastq -- see entry points below
    path: /abs/path/to/data
    busco_lineage: fungi_odb12   # optional, default shown
    organelle: none               # optional, GetOrganelle db or 'none'
    subsample: {coverage: 100, genome_size: 40mb}  # optional coverage cap
    porechop: true                 # optional, adapter trimming (default on)
    filter:                        # optional, PER-SAMPLE grid override
      min_q:   [20]
      min_len: [1000]

dorado:
  version: "2.1.2"          # or: path: /my/own/dorado/bin/dorado

dorado_correct:
  enabled: false             # hifiasm --ont does its own correction; see below
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

### The read-filtering grid, and the per-sample override
`chopper` filters reads at every `(min_q, min_len)` combination in `filter:`,
and each cell gets assembled independently with `hifiasm`. **This top-level
`filter:` is the default for every sample that doesn't set its own** (see
`sample_a` in `config/config.yaml`'s example, which overrides it to a single
cell). If you give one sample a much larger or smaller grid than another,
give it its own `filter:` block -- otherwise editing the shared default for
one sample silently changes what every *other* sample's `contig_count`/
`highest_busco`/`assembly_stats` requires too, including grid cells that
sample never actually ran. This is a real trap, not a hypothetical one: it's
exactly what happened during this pipeline's own development, when one
sample's grid was widened for a separate experiment and silently changed
what another, already-finished sample's `contig_count` demanded next.

**Two independent rules pick a winner across that grid, in parallel:**
- `contig_count` picks the fewest-contigs assembly. For a repeat-rich
  genome, don't take this at face value: check the winner's total length in
  `assembly_stats.tsv` against your expected genome size -- fewer contigs
  can mean collapsed repeats rather than a genuinely better assembly.
- `highest_busco` runs BUSCO on *every* grid cell's assembly and picks the
  one with the highest completeness ("C") percentage instead -- an
  orthogonal, biology-grounded signal that doesn't share `contig_count`'s
  collapsed-repeat blind spot. **Cost warning:** this runs BUSCO once per
  grid cell, not once total -- a 3x4=12-cell grid means 12 BUSCO runs just
  for this selection, on top of the two more (one per selector) that run
  again afterwards on each winner's final, possibly-polished assembly.

Both winners go through the rest of the pipeline independently -- polishing
if a BAM is available, a final BUSCO run, stats, and a deliverable -- so
every sample ends up with two directly comparable results:
`results/<sample>/lowest_contig/` and `results/<sample>/highest_busco/`.
Neither is automatically "more correct"; compare their `summary.txt` files
(contig count, N50, total length vs. expected genome size, and BUSCO score)
and use your judgement.

### `dorado_correct` vs hifiasm's own correction
`hifiasm --ont` does its own ONT-specific read correction, and that's the
default (`dorado_correct.enabled: false`). Only turn `dorado_correct` on if
you specifically want to experiment with dorado's correction instead: on
one large test sample (N50 ~5.2 kb) a real GPU run of `dorado correct`
discarded ~85% of reads (7.27 Gbp in -> 1.11 Gbp out), dropping ~30x raw
coverage to ~4.6x effective coverage -- far too little for hifiasm to
assemble anything useful from. `dorado correct` may still be worth using on
a dataset with a longer read N50, where it discards proportionally less.

### BUSCO lineage
Set `busco_lineage` per sample to match the organism -- e.g. a
`stramenopiles_odb*` dataset for an oomycete like *Phytophthora*, which is
not a fungus despite this pipeline's name. `busco --list-datasets` (inside
the BUSCO conda env) lists what's available. Scores aren't comparable across
different dataset versions (odb10 vs odb12), so don't mix them when
comparing runs.

## Running the pipeline

**First run only -- build the conda environments** (Snakemake itself needs
almost no resources, but building ~10 tool environments does):
```
snakemake --profile profile --conda-create-envs-only
```

**Then run the workflow.** `slurm_submit.sbatch` wraps this in a single
lightweight SLURM job (1 CPU, 1G -- each *rule* still gets its own
appropriately-sized job via the executor plugin):
```
sbatch slurm_submit.sbatch
squeue --me      # watch it
```
Alternatively, since the orchestrator itself is cheap and just submits/polls,
run it directly on the login node inside `tmux` instead of wrapping it in an
sbatch job -- this is the more common way to drive a long multi-day run
interactively:
```
tmux new -s my_run
conda activate snakemake_run
snakemake --profile profile
# detach with Ctrl-b d; reattach later with: tmux attach -t my_run
```

**Multiple samples configured at once:** `rule all` builds every sample in
`config.yaml` in full. While iterating on one sample, target its outputs
explicitly instead of running bare `snakemake --profile profile`, e.g.:
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
as soon as their consumers finish, which also deletes the evidence you'd
want to inspect. Re-run with `--notemp` to keep everything while debugging.

## Output

```
results/<sample>/
  read_stats.tsv                      # seqkit stats for the raw (and trimmed) reads -- shared, sample-level

  lowest_contig/                      # winner selected by fewest contigs
    summary.txt                       #   human-readable: reads, grid, coverage, BUSCO
    assembly_stats.tsv                #   seqkit stats for every grid candidate + this winner
    <sample>_lowest_contig_final.fasta  # THE deliverable for this selector -- polished if a BAM was available

  highest_busco/                      # winner selected by highest BUSCO completeness
    summary.txt
    assembly_stats.tsv
    <sample>_highest_busco_final.fasta  # THE deliverable for this selector

logs/                   # per-rule, per-sample logs
snake_log/               # the Snakemake orchestrator's own run logs
```

Both selectors run all the way through independently (see "The
read-filtering grid" above) -- there is no single "the" final assembly,
there are two: compare `lowest_contig/summary.txt` and
`highest_busco/summary.txt` and pick whichever actually looks better for
your genome.

`data/` holds intermediates. Most of it is `temp()` and gets cleaned up
automatically once nothing downstream needs it any more -- don't expect
`data/contig/`, `data/dorado_polish/`, `data/busco_grid/*/BUSCO` (BUSCO on
every grid cell, used only for the highest_busco selection), or
`data/busco/*/BUSCO`'s full tables to persist after a run finishes; that
data is either duplicated in `results/` already or (for BUSCO's full
tables/logs) considered disposable once the short summary is captured in
`summary.txt`. The two things that persist deliberately:
`data/dorado_basecall/<sample>.bam` (`protected()` -- basecalling is
expensive to redo) and, if a sample sets `organelle:`,
`data/getorganelle/<sample>/Mitochondria` (there's no `results/` copy of
that path yet -- it's the only copy of the organelle assembly, so it's
*not* `temp()`).

## Repo layout
- `Snakefile` -- entry-point handling, input resolvers, adaptive resource
  helpers, and the rule includes/targets. Commented throughout; read it
  before adding a new rule or entry point.
- `workflow/rules/*.smk` -- one file per pipeline stage, numbered in
  roughly the order they run in the DAG. Each has a header comment
  explaining what it does and when it's reached.
- `workflow/envs/*.yml` -- one conda env per tool.
- `workflow/scripts/summarize.py` -- builds
  `results/<sample>/<selector>/summary.txt`.
- `workflow/scripts/pick_highest_busco.py` -- picks the `highest_busco`
  selector's grid winner (see `rule highest_busco`).
- `config/config.yaml`, `profile/config.yaml` -- see above.
- `scripts/` -- separate, standalone scripts used to make the paper's
  figures (phylogeny, circos plots, telomere/mitochondria checks, ...).
  Barebones and not part of the Snakemake DAG; see `scripts/README.md`.

## Cluster notes (AAU BioCloud)
Learned from running this pipeline on BioCloud; useful if you're adapting
`profile/config.yaml`/`slurm_submit.sbatch` to a different cluster too.

- **`--partition` has no effect here.** BioCloud assigns partitions
  automatically from each job's memory-per-CPU ratio and node features --
  that's why no rule or profile setting here specifies one.
- **GPU:** exactly one node, `bio-node10` -- 1x NVIDIA A10 (Ampere, CC 8.6,
  23 GB VRAM), 64 threads, 256 GB RAM. Request it with `--gres=gpu:a10:1`
  (BioCloud's GPU_GRES in the Snakefile), not via `--partition`.
- **Interactive jobs land on the `interactive` partition** (`bio-node11`, no
  GPU) unless you explicitly pass `--gres` for GPU work interactively.
- **Fat CPU nodes** (`zen5x` ~2.3 TB RAM/288 threads, `zen3x` ~2 TB) are
  better than the GPU node for memory-hungry non-GPU work like hifiasm or
  dorado correct's overlap stage.
- **Compute nodes have internet access** -- conda env creation and dorado's
  download both work from inside a job, not just the login node.
- **The SLURM account auto-guess can fail** (`sacct: invalid option -- '1'`)
  without breaking anything. If jobs start getting rejected for account
  reasons, set an explicit `slurm_account` under `default-resources` in
  `profile/config.yaml`.
- **Run the Snakemake orchestrator on the login node, not inside a SLURM
  job.** Snakemake warns about this; the orchestrator itself barely uses
  any resources (it only submits jobs and polls), so there's no benefit to
  wrapping it in its own allocation, and doing so can make its own SLURM
  job look like it's hung.

## Known gotchas
- **`temp()` deletes your debugging evidence.** Re-run with
  `snakemake --profile profile --notemp` while debugging a failing rule so
  its intermediates stick around for inspection.
- **A missing declared input crashes DAG *construction*, even for an
  already-complete, up-to-date target.** Snakemake resolves every reachable
  input's producer rule while building the DAG, before it checks whether
  the consuming output actually needs rebuilding. If you delete a sample's
  raw data (or another intermediate) once its `results/` are final, a later
  plain `snakemake --profile profile` run touching that sample can crash
  with `MissingInputException` trying to plan how to regenerate the deleted
  file -- regardless of `rerun-triggers` and regardless of whether the
  target you actually asked for is already fine. Fix: leave a 0-byte
  placeholder at the deleted path with an old timestamp
  (`touch -d "<old date>" <deleted path>`) -- that satisfies DAG
  construction without needing the real data back.
- **`--until <rule>`** stops the pipeline at a chosen rule, useful for
  building just far enough to debug one stage.
- **Changing a conda env `.yml` creates a new, separately-hashed env** the
  next time that rule runs. Old ones accumulate under `.snakemake/conda/`;
  `rm -rf .snakemake/conda` after a round of version bumps to avoid keeping
  several stale environments around.
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

## Notes
- If you only have a basecalled BAM (no raw pod5), set `type: bam` --
  basecalling is skipped but polishing still runs.
- If you only have FASTQ (no BAM at all), set `type: fastq` -- both
  basecalling and polishing are skipped, since polishing needs a BAM's move
  table. The pipeline still filters, assembles, and reports on the winning
  grid cell.
- Do NOT set your conda channel priority to `strict`, or you'll end up in
  dependency hell resolving these environments.
- `.gitignore` excludes everything generated at runtime (`data/`, `results/`,
  `logs/`, `snake_log/`, `resources/`, efficiency reports, ...). Check
  `git status` before `git add -A` -- it's easy to accidentally stage
  generated junk.

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
