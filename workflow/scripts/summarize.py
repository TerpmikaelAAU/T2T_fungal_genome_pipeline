"""Build a single human-readable summary per sample.

Pulls together seqkit read stats, seqkit assembly stats, the BUSCO short
summary, and a coverage estimate derived from corrected bases / assembly size.
Written to be forgiving: if BUSCO output cannot be located the rest of the
report is still produced, rather than failing the whole run at the last step.
"""

import glob
import os
import re

snakemake = snakemake  # noqa: F821  (injected by Snakemake)


def read_tsv(path):
    """seqkit stats -T output -> list of dicts."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def as_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def label(path):
    """Turn a data/ path into something readable in the report."""
    base = os.path.basename(path)
    parent = os.path.basename(os.path.dirname(path))
    if "porechopped" in path:
        return "after adapter trimming"
    if "hifiasm" in path:
        return f"assembly {parent}"
    if "contig" in path:
        return "selected assembly"
    if "decompressed" in path or "samtools/Fastq" in path:
        return "raw reads"
    # Whatever the user pointed `samples: <name>: path:` at directly (a
    # non-gz fastq, or a fastq.gz small enough that read_stage_files skipped
    # the decompress-once rule for something else's sake -- see Snakefile
    # get_raw_fastq()).
    return "raw reads"


def find_busco(busco_dir):
    """BUSCO's short summary lives at an unpredictable depth; search for it."""
    hits = glob.glob(os.path.join(busco_dir, "**", "short_summary*.txt"),
                     recursive=True)
    return hits[0] if hits else None


reads = read_tsv(snakemake.input.reads)
asm = read_tsv(snakemake.input.asm)

# --- raw bases, for the coverage estimate ----------------------------------
# Filtering/correction now fork into a (min_q, min_len) grid (see
# config.yaml `filter:`), so there is no longer a single corrected-reads
# file to measure against -- this is coverage of the UNFILTERED input,
# a much looser upper bound than the old corrected-bases estimate.
raw_bases = 0
for row in reads:
    f = row.get("file", "")
    if "porechopped" not in f and "hifiasm" not in f and "contig" not in f:
        raw_bases = as_int(row.get("sum_len"))
        break

# --- the selected assembly ------------------------------------------------
best = None
for row in asm:
    if "contig/" in row.get("file", ""):
        best = row
best_len = as_int(best.get("sum_len")) if best else 0

coverage = (raw_bases / best_len) if best_len else 0.0

out = []
out.append("=" * 70)
out.append(f"  SUMMARY: {snakemake.params.sample}")
out.append("=" * 70)
out.append("")
out.append(f"Polished          : {'yes' if snakemake.params.polished else 'no (FASTQ entry point)'}")
out.append(f"Final assembly    : {snakemake.input.final}")
out.append("")

out.append("-" * 70)
out.append("READS")
out.append("-" * 70)
out.append(f"{'stage':<32}{'reads':>12}{'bases':>16}{'N50':>12}")
for row in reads:
    out.append(f"{label(row.get('file','')):<32}"
               f"{as_int(row.get('num_seqs')):>12,}"
               f"{as_int(row.get('sum_len')):>16,}"
               f"{as_int(row.get('N50')):>12,}")
out.append("")

n_candidates = sum(1 for row in asm if "contig/" not in row.get("file", ""))

out.append("-" * 70)
out.append("ASSEMBLIES  (one per min_q/min_len grid cell; lowest contig count is selected)")
out.append("-" * 70)
if n_candidates == 1:
    out.append("Note: only 1 grid cell was configured for this sample -- selected")
    out.append("assembly is that cell's output, not a winner across a grid.")
    out.append("")
out.append(f"{'candidate':<32}{'contigs':>12}{'total bp':>16}{'N50':>12}")
for row in asm:
    out.append(f"{label(row.get('file','')):<32}"
               f"{as_int(row.get('num_seqs')):>12,}"
               f"{as_int(row.get('sum_len')):>16,}"
               f"{as_int(row.get('N50')):>12,}")
out.append("")

out.append("-" * 70)
out.append("COVERAGE")
out.append("-" * 70)
out.append(f"Raw bases         : {raw_bases:,}")
out.append(f"Assembly size     : {best_len:,}")
out.append(f"Estimated coverage: {coverage:.1f}x")
out.append("  (raw bases / assembly size -- an upper bound, not a mapped depth)")
out.append("")

out.append("-" * 70)
out.append(f"BUSCO  ({snakemake.params.lineage})")
out.append("-" * 70)
summary_file = find_busco(str(snakemake.input.busco))
if summary_file:
    with open(summary_file) as fh:
        for line in fh:
            if re.search(r"C:|Complete|Fragmented|Missing|Total BUSCO", line):
                out.append("  " + line.strip())
else:
    out.append("  short_summary*.txt not found under "
               f"{snakemake.input.busco}")
    out.append("  (check the BUSCO rule's -o path; see the run log)")
out.append("")

with open(snakemake.output.a, "w") as fh:
    fh.write("\n".join(out) + "\n")