"""Pick the grid cell with the fewest contigs -- the `lowest_contig` selector
(6_00_lowest_contig_count.smk). Only cells of at least the sample's
min_assembly_mb count (see assembly_size.py); ties go to the first path in
sort order.
"""

import shutil
import sys

snakemake = snakemake  # noqa: F821  (injected by Snakemake)
sys.path.insert(0, snakemake.scriptdir)
from assembly_size import eligible  # noqa: E402

fasta_files = list(snakemake.input.fasta_files)
with open(snakemake.log[0], "w") as log:
    candidates = eligible(fasta_files, snakemake.params.min_bp, log)
    if not candidates:
        raise ValueError("No valid contig counts found -- every grid point "
                         "produced an empty or missing assembly.")
    best = min(candidates, key=lambda c: (c[1], c[0]))[0]
    if len(fasta_files) == 1:
        print("Note: only 1 grid cell configured for this sample -- picked it "
              "directly, nothing to compare against.", file=log)
    print(f"picked: {best}", file=log)

shutil.copy(best, str(snakemake.output.a))
