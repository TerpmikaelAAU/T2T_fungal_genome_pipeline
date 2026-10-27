# T2T fungal genome assembly pipeline

> **This is a work in progress.**

Paper: "Gapless near Telomer-to-Telomer Assembly of Neurospora intermedia,
Aspergillus oryzae, and Trichoderma asperellum from Nanopore Simplex Reads",
https://doi.org/10.3390/jof11100701. The `v1.0.0` tag is the pipeline the
paper describes. Later versions add flexible entry points, a read-filtering
grid and resources that scale with input size, and have also been run on a
~240 Mb genome that is ~75% transposable elements.

## What it does
From Nanopore reads (pod5, basecalled BAM or FASTQ, or a public SRA/ENA run):

1. Downloads (sra), basecalls (pod5) and optionally trims adapters.
2. Filters the reads over a grid of quality/length cutoffs and assembles
   each cell with `hifiasm --ont`.
3. Picks winners across the grid in independent ways: fewest contigs,
   highest BUSCO, and (with telomeres on) most telomere-to-telomere contigs.
4. Polishes each winner with `dorado polish` when the reads allow it.
5. Reports BUSCO (before and after polishing), read/assembly stats,
   telomeres, a Tiara contamination screen and an HTML report per sample.
6. Optional: organelle genome (Flye + GetOrganelle), annotation the BAGS
   way, and a UFCG phylogeny with one NCBI genome per fungal genus.

Each rule instance is its own SLURM job, submitted by Snakemake.

## Requirements
- Nanopore reads, ideally R10.4 chemistry or newer, at good coverage (the
  paper used 100x+ at Q10).
- AAU BioCloud (https://cmc-aau.github.io/biocloud-docs/), or another
  cluster with Snakemake 9 and the SLURM executor plugin.

## Installation
```
conda env create -n snakemake_run -f Snakemake_env.yml
conda activate snakemake_run
```
This pins Snakemake 9.26.1 and the SLURM executor plugin; Snakemake 7/8
won't work. Tools run in containers that are pulled on first use, and dorado
is downloaded on first use, so there is nothing else to install.
`profile/config.yaml` is set up for BioCloud.

`config/config.yaml` ships pointing at `example_data/`: tiny fake reads for
trying the scheduling and containers before a real run (see
`example_data/README.md`). They are too few to assemble, so the run stops
at the selectors (`contig_count`, `highest_busco`, `most_t2t`) with "No valid
..." / "No non-empty assembly" errors.

## Configuring a run
All settings are in `config/config.yaml`. Each entry under `samples:` is one
genome and names its `results/` folder:

```yaml
samples:
  my_fungus:
    type: bam
    path: /absolute/path/to/my_fungus.bam
    porechop: false
  rhizopus_arrhizus:
    type: sra
    accession: SRR36829119
```

Keep cluster paths out of commits: CI dry-runs the committed config and
fails on paths that only exist on the cluster.

### Per-sample keys
| Key | Default | What it does |
|---|---|---|
| `type` | required | `pod5`, `bam`, `fastq` or `sra`; see "Entry points". |
| `path` | required, except for `sra` | A directory of `.pod5` files, ONE basecalled `.bam` (merge several with `samtools merge`), or one `.fastq`/`.fastq.gz`. Absolute, or relative to the repo root. |
| `accession` | required for `sra` | A run accession (`SRR`/`ERR`/`DRR`), or a list when one sample spans several runs. |
| `busco_lineage` | `fungi_odb12` | BUSCO dataset. Match it to the organism (e.g. `stramenopiles_odb12` for an oomycete). Don't compare scores across odb versions. |
| `min_assembly_mb` | global | Smallest grid assembly worth picking. |
| `porechop` | `true` | porechop_abi adapter trimming. Slow on big BAMs, and dorado already trims while basecalling. |
| `filter` | global `filter` | This sample's own grid; see "The grid". |
| `ultralong` | global `ultralong.enabled` | Also give hifiasm the reads above `ultralong.min_len` as `--ul`. |
| `dorado_correct` | global `dorado_correct.enabled` | dorado read correction instead of hifiasm's own; see "The grid". |
| `subsample` | off | rasusa coverage cap, e.g. `{coverage: 100, genome_size: 40mb}`. Using all reads is recommended. |
| `organelle` | `false` | Also recover an organelle genome with Flye + GetOrganelle. |
| `organelle_type` | required if `organelle` | GetOrganelle's `-F` type; fungi want `fungus_mt`. |

### Global settings
| Setting | Shipped | What it does |
|---|---|---|
| `dorado.version` | `"2.1.2"` | dorado release fetched on first use. Use `path:` instead to point at your own install. |
| `dorado_polish.device` | `"cuda:all"` | `"cpu"` polishes on a normal node when the GPU node is down or busy. Much slower. |
| `dorado_correct.enabled` | `false` | Default for each sample's `dorado_correct`. |
| `dorado_correct.min_q`, `min_len` | `10`, `10000` | The one cutoff applied to reads going into dorado correct. |
| `dorado_correct.index_size` | `"4G"` | Smaller means less RAM per block but more blocks. |
| `dorado_correct.max_blocks` | `200` (not in the shipped config) | Stop if dorado reports more blocks than this, e.g. after a typo like `index_size: "4"`. |
| `dorado_correct.inference_device` | `"cuda:all"` | `"cpu"` keeps inference off the GPU (overlap is always CPU). |
| `ultralong.enabled` | `true` | Default for each sample's `ultralong`. |
| `ultralong.min_q`, `min_len` | `10`, `50000` | Cutoff for the `--ul` reads. Check there is coverage above it: one large dataset had only ~3.5x above 20 kb. |
| `hifiasm_min_input_mb` | `1` | Cells with less read data skip hifiasm and get an empty placeholder, because hifiasm crashes with an illegal instruction on a near-empty FASTQ (and retrying with more memory only repeats the crash). |
| `min_assembly_mb` | `10` | Assemblies smaller than this are never picked and get no BUSCO run (BUSCO can crash on a gene-less assembly). Keep it well below the expected genome size; `0` turns it off. |
| `prefilter.enabled`, `min_q`, `min_len` | `true`, `10`, `10000` | Filter the full read set once and let every stricter cutoff read that smaller file. Same reads per cell either way; best set to the grid's loosest cell. |
| `filter.min_q`, `min_len` | `[10, 15, 20]`, `[10000, 15000, 20000, 25000, 30000]` | The default grid. |
| `telomere.enabled` | `true` | Telomere motif, report and the `most_t2t` selector; see "Telomeres". Off if the `telomere:` block is missing. |
| `telomere.top_motifs` | `10` | Candidate motifs from `tidk explore` tried with `tidk search`. |
| `telomere.explore_min_len`, `explore_max_len` | `5`, `12` | Repeat unit lengths (bp) `tidk explore` tries. |
| `telomere.explore_threshold` | `10` | Units in a row that count as a repeat run. tidk's default of 100 misses short fungal telomeres. |
| `telomere.window`, `min_repeats` | `2000`, `10` | An end is telomeric with at least `min_repeats` motif copies in its last `window` bp. |
| `telomere.rerun_hifiasm` | `true` | Re-assemble every cell with `hifiasm --telo-m <motif>` when the motif is clear (below). Doubles the hifiasm cost. |
| `telomere.rerun_min_fraction`, `rerun_min_ends` | `0.8`, `4` | The motif is clear when it covers at least this fraction of all telomeric ends found and at least this many ends. |
| `telomere.fallback_motif` | `TTAGGG` | Used for the counts when no candidate finds a telomere; hifiasm isn't re-run then. |
| `contamination.enabled`, `min_len` | `true`, `3000` | Tiara screen; contigs shorter than `min_len` aren't classified. Off if the `contamination:` block is missing. |
| `report.enabled` | `true` | Per-sample HTML report, new runs only; see "The report". |
| `phylogeny.enabled` | `false` | Genus-level tree; see "Phylogeny". |
| `phylogeny.taxon` | `"Fungi"` | Any NCBI taxon name or taxid, e.g. `Ascomycota`. Changing it means deleting `data/phylogeny/`. |
| `bags.enabled` | `false` | Annotation; see "Annotation". |
| `bags.geneml_version` | `"1.1.0"` | geneML version installed from PyPI. |

### Entry points
| `type:` | Runs |
|---|---|
| `pod5` | basecall → filter → assemble → polish |
| `bam` | filter → assemble → polish |
| `fastq` | filter → assemble → polish only if the reads carry dorado's read-group tags |
| `sra` | download → filter → assemble, never polished |

**Polishing** needs to know which basecaller model made the reads. A dorado
BAM records it in its header. A FASTQ only has it if every header line
carries a tag like `RG:Z:<run id>_dna_r10.4.1_e8.2_400bps_sup@v5.0.0`; check
with `zcat -f reads.fastq.gz | head -1`. SRA rewrites read headers and drops
the tags. The pipeline prints at startup whether each FASTQ sample will be
polished. Reads with a move table (`--emit-moves`, always on for pod5) get
dorado's more accurate polishing model.

**sra samples** are downloaded from ENA as `data/sra/<sample>.fastq.gz`.
- ENA is slow (1-2 MB/s, a few hours for 10 GB) and drops long transfers.
  The download resumes from `data/sra/<sample>.fastq.gz.partial.<file>`, so
  a rerun carries on where it stopped.
- Each file is md5-checked and downloaded again once on a mismatch. If it
  still fails, download the run with `prefetch` + `fasterq-dump` and use it
  as a `fastq` sample.
- `data/sra/` is kept (about 1 GB per Gbase); delete it once `results/` are
  final. Check the run at https://www.ebi.ac.uk/ena/browser/view/<accession>
  first: non-Nanopore runs are refused.

**Duplex data** from `dorado duplex` (bam or fastq) needs no config. The
simplex parents of duplex reads are dropped before assembly, so each
molecule counts once (by `dx:i:-1` in a BAM, by the `<id1>;<id2>` read name
in a FASTQ, so the names must be dorado's own). `dorado polish` rejects
duplex reads, so polishing gets all simplex reads, the parents included.

### The grid
`chopper` filters the reads at every `(min_q, min_len)` in `filter:` and
each cell is assembled on its own. Give a sample its own `filter:` instead
of editing the global one: changing the global grid changes what every other
sample's selectors expect, including cells they never ran.

The selectors:
- `lowest_contig` (fewest contigs). For a repeat-rich genome, fewer contigs
  can mean collapsed repeats; check the total length against the expected
  genome size.
- `highest_busco` runs BUSCO on every cell, one run per cell on top of the
  final BUSCO runs.
- `most_t2t` (with telomeres on); see "Telomeres".

Every winner is BUSCO'd and reported on its own under
`results/<sample>/<selector>/` (an assembly picked by several selectors is
polished once). None is automatically better; compare
their `summary.txt`. If no cell reaches `min_assembly_mb`, the selectors fall
back to every non-empty cell and say so in their logs (`logs/select/`,
`logs/telomere/` for `most_t2t`).

**`dorado_correct`** is off by default because `hifiasm --ont` corrects ONT
reads itself, and on a ~5 kb N50 dataset dorado correct discarded ~85% of
the reads. When on, it runs once per sample at the `dorado_correct` cutoff,
and only the grid's `min_len` values apply afterwards. Its output has no
qualities, so it is padded with Q30: hifiasm drops reads below Q10.

## Telomeres
On in the shipped config; the keys are in "Global settings".
1. `tidk explore` lists candidate repeat units near contig ends of the
   fewest-contig first-pass assembly; each of the top `top_motifs` is tried
   with `tidk search`, and the one making most ends telomeric wins
   (`results/<sample>/telomere_motif_candidates.tsv`).
2. If the motif is clear (`rerun_min_fraction`, `rerun_min_ends`), every
   cell is re-assembled with `hifiasm --telo-m <motif>` and the selectors
   choose across both runs.
3. `most_t2t` picks the cell with the most telomere-to-telomere contigs.
4. Every final assembly gets `telomeres.tsv`, a tidk plot and a TELOMERES
   section in `summary.txt`.

## Contamination screen
Tiara labels each contig of every final assembly as eukarya, bacteria,
archaea, organelle or unknown; `summary.txt` warns about bacterial/archaeal
contigs. It only reports, nothing is removed. Tiara rather than NCBI FCS-GX
because FCS-GX wants a ~500 GB database in RAM, and NCBI runs it on
submission anyway. Tiara isn't on bioconda and needs Python ≤ 3.9, so the
first run installs it into `resources/tiara-<version>/` inside
`python:3.9-slim`.

## Annotation
`bags.enabled: true` runs the [BAGS](https://github.com/TerpmikaelAAU/BAGS)
steps on every final assembly: geneML genes, antiSMASH on those genes, and
BUSCO in protein mode (how complete the annotation is). One row per
assembly goes in `results/annotation_summary.tsv`. geneML is only on PyPI,
so it is installed once into `resources/geneml-<version>/`; the antiSMASH
databases are downloaded once into `resources/antismash_databases/`.

## Phylogeny
`phylogeny.enabled: true` builds one [UFCG](https://ufcg.steineggerlab.com)
core-gene tree (FastTree) over every sample's final assemblies plus one
NCBI genome per genus under `taxon` (the NCBI reference genome, else the
best assembly level, else the highest contig N50).
- All of Fungi is a few thousand genomes, ~100 GB, and one UFCG job per
  genome. `data/phylogeny/` is kept so this happens once.
- The reference jobs start straight away, rank below everything else and use
  at most `ufcg_reference_slots` (10 in `profile/config.yaml`) of the 30
  job slots, so they never hold up the assemblies.
- Genomes UFCG can't profile are left out with a warning in
  `logs/phylogeny/`.
- `export NCBI_API_KEY=<key>` raises NCBI's rate limit.
- Genus/family level only, and the markers are fungal, so not for oomycetes.
- Skipped when every sample is from `example_data/`.

## Containers
Each rule's `container:` names a pinned biocontainers image, pulled on first
use. Exceptions: NCBI `datasets` (StaPH-B's image; it left bioconda), UFCG
(the official image, which ships the core-gene database), and geneML and
Tiara (installed into plain Python images, above). `workflow/envs/*.yml` are
reference only.

`scripts/` holds the standalone scripts behind the paper's figures; they are
not part of the workflow.

## Running
On the login node inside `tmux`:
```
tmux new -s my_run
conda activate snakemake_run
snakemake --profile profile
```
or as one small SLURM job with `sbatch slurm_submit.sbatch` (edit the conda
env name and `--mail-user` in it first; see "Cluster notes" for the
drawback).
- Dry-run first with `-n`, especially after pulling changes into a running
  project: most intermediates are `temp()`, and a new consumer of one can
  rebuild it, which can mean re-running hifiasm for every cell.
- Build one sample only by naming its targets, e.g.
  `results/my_sample/lowest_contig/summary.txt`.
- After a crash: `snakemake --profile profile --unlock`.
- Debugging: `--notemp` keeps the intermediates that `temp()` would delete.

## Output
```
results/<sample>/
  report.html                       # charts of reads, grid, k-mer spectra, BUSCO and picks
  read_stats.tsv                    # reads at each stage and per grid cell
  polish_groups.tsv                 # polished samples: selectors that picked the same assembly
  telomere_motif_candidates.tsv
  <selector>/                       # lowest_contig, highest_busco, most_t2t
    summary.txt                     # reads, grid, coverage, BUSCO, telomeres, contamination
    assembly_stats.tsv
    <sample>_<selector>_final.fasta # the deliverable
    telomeres.tsv, telomeres/
    contamination/tiara.tsv
    annotation/                     # bags only
results/annotation_summary.tsv      # bags only
results/phylogeny/                  # genus_tree.nwk, alignment, tip_labels.tsv
logs/                               # per rule and sample; logs/hifiasm/ has hifiasm's own output
```

Most of `data/` is deleted once nothing needs it. Kept on purpose:
`data/dorado_basecall/` (basecalling is expensive), `data/getorganelle/`
(the only copy), `data/sra/`, `data/phylogeny/`, and `data/report/` (about
10 bytes per read) so the report can be redrawn.

### The report
`results/<sample>/report.html` is one self-contained page: read length and
quality, coverage per cutoff, hifiasm's k-mer spectrum per cell, BUSCO, and
a scoreboard of all cells with each selector's pick.

It is built for **new runs only**. It needs the grid's input reads and
BUSCO results, which are `temp()`; asking for it on a finished run would
rebuild them, hifiasm included. So a sample gets a report when it hasn't
started, or while those files still exist or are still to come; an
existing `report.html` always stays a target. Cells
assembled before `logs/hifiasm/` existed show "no histogram".

## Cluster notes (AAU BioCloud)
- `--partition` does nothing: BioCloud picks the partition from the
  memory-per-CPU ratio.
- One GPU node, `bio-node10` (1x A10, 23 GB VRAM, 64 threads, 256 GB RAM);
  request it with `--gres=gpu:a10:1`.
- Interactive jobs land on `bio-node11` (no GPU) unless you pass `--gres`.
- The fat CPU nodes (`zen5x` ~2.3 TB, `zen3x` ~2 TB) suit hifiasm and dorado
  correct's overlap step better than the GPU node.
- Compute nodes have internet access.
- The SLURM account guess fails (`sacct: invalid option -- '1'`) harmlessly.
  If jobs get rejected over the account, set `slurm_account` under
  `default-resources` in `profile/config.yaml`.
- Snakemake warns when it runs inside a SLURM job (`slurm_submit.sbatch`).
  It only submits and polls, so the login node is enough, and inside a job it
  can look hung.

## Known gotchas
- Deleting a sample's raw data after `results/` are final makes later runs
  crash with `MissingInputException`. Leave an old 0-byte file in its place:
  `touch -d "2020-01-01" <deleted path>`.
- Changing a rule's `container:` pulls a new image; old ones pile up in
  `.snakemake/singularity/`.
- If you do use `workflow/envs/*.yml`, don't set conda channel priority to
  `strict`; they won't resolve.

## How to cite
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

When using a version after `v1.0.0`, please also give the tag you ran and
link this repository.
