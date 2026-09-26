"""Pick the grid cell with the most telomere-to-telomere contigs -- the
`most_t2t` selector (6_02_most_t2t.smk), alongside fewest contigs and
highest BUSCO.

Ties go to more telomeric ends, then fewer contigs. Empty placeholder
cells are skipped; if no cell has any telomere at all, this still picks
by those tie-breakers (i.e. fewest contigs) rather than failing.
"""

import os
import shutil
import sys

snakemake = snakemake  # noqa: F821  (injected by Snakemake)
sys.path.insert(0, snakemake.scriptdir)
from telomeres import contig_ends, read_search, tally  # noqa: E402

searches = list(snakemake.input.searches)
fasta_files = list(snakemake.input.fasta_files)
if len(searches) != len(fasta_files):
    raise ValueError("searches and fasta_files are out of sync -- both must "
                     "follow grid_assemblies() order.")

best_key, best_fasta = None, None
with open(snakemake.log[0], "w") as log:
    print("cell\tt2t_contigs\ttelomeric_ends\tcontigs", file=log)
    for search, fasta in zip(searches, fasta_files):
        if not os.path.getsize(fasta):
            continue
        ends, t2t, n = tally(contig_ends(read_search(search), snakemake.params.min_repeats))
        if n == 0:
            continue
        print(f"{fasta}\t{t2t}\t{ends}\t{n}", file=log)
        key = (t2t, ends, -n)
        if best_key is None or key > best_key:
            best_key, best_fasta = key, fasta
    print(f"picked: {best_fasta}", file=log)

if best_fasta is None:
    raise ValueError("No non-empty assembly across the grid -- nothing to pick.")
shutil.copy(best_fasta, str(snakemake.output.a))
