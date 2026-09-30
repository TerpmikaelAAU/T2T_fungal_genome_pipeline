"""Pick the assembly with the highest BUSCO completeness ('C' percentage)
across a sample's read-filtering grid -- an alternative to contig_count's
fewest-contigs pick (6_00_lowest_contig_count.smk), using an independent,
biology-grounded signal instead of contig count alone.

Written to be forgiving like contig_count: an empty or missing candidate
scores as absent, not a crash, and only fails if EVERY candidate is unusable.
Cells smaller than the sample's min_assembly_mb don't count (see
assembly_size.py).
"""

import glob
import os
import re
import shutil
import sys

snakemake = snakemake  # noqa: F821  (injected by Snakemake)
sys.path.insert(0, snakemake.scriptdir)
from assembly_size import eligible  # noqa: E402


def complete_pct(busco_dir):
    """Parse the 'C:97.7%[...]' figure out of a BUSCO short_summary file.

    Returns None if no summary can be found or parsed (treated the same as
    a missing/empty candidate -- skipped, not fatal)."""
    hits = glob.glob(os.path.join(busco_dir, "**", "short_summary*.txt"),
                      recursive=True)
    if not hits:
        return None
    with open(hits[0]) as fh:
        match = re.search(r"C:([0-9.]+)%", fh.read())
    return float(match.group(1)) if match else None


busco_dirs = list(snakemake.input.busco_dirs)
fasta_files = list(snakemake.input.fasta_files)

if len(busco_dirs) != len(fasta_files):
    raise ValueError(
        f"busco_dirs ({len(busco_dirs)}) and fasta_files ({len(fasta_files)}) "
        "are out of sync -- both must expand() over the same grid, in the "
        "same order, for pairing by position to be valid."
    )

busco_of = dict(zip(fasta_files, busco_dirs))
best_score, best_fasta = None, None
with open(snakemake.log[0], "w") as log:
    for fasta, _, _ in eligible(fasta_files, snakemake.params.min_bp, log):
        score = complete_pct(busco_of[fasta])
        print(f"{fasta}\tBUSCO C: {score}", file=log)
        if score is None:
            continue
        if best_score is None or score > best_score:
            best_score, best_fasta = score, fasta
    print(f"picked: {best_fasta}", file=log)

if best_fasta is None:
    raise ValueError(
        "No valid BUSCO score found across the grid -- every grid point "
        "produced an empty assembly or an unparseable BUSCO summary."
    )

if len(fasta_files) == 1:
    print("Note: only 1 grid cell configured for this sample -- picked it "
          "directly, nothing to compare against.", file=sys.stderr)

shutil.copy(best_fasta, str(snakemake.output.a))
