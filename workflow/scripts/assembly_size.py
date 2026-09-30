"""Shared by the grid selectors (pick_lowest_contig.py, pick_highest_busco.py,
pick_most_t2t.py): which grid cells are big enough to be picked at all.

A cell whose total assembly size is below the sample's min_assembly_mb (see
config.yaml) is a fragment, not a genome -- on sparse data a single 50 kb
contig would otherwise "win" fewest contigs. If NO cell reaches the minimum,
every non-empty cell stays eligible (with a warning) so the run still
finishes and the numbers can be inspected.
"""

import os


def fasta_stats(path):
    """(number of contigs, total bp) of a FASTA file."""
    n = bp = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
            else:
                bp += len(line.strip())
    return n, bp


def eligible(fastas, min_bp, log):
    """Non-empty candidates at least min_bp long, as (fasta, contigs, bp)
    tuples in input order; all non-empty ones if none is that long. Writes
    every candidate's size and the outcome to `log` (an open file)."""
    non_empty = []
    print("cell\tcontigs\tbp", file=log)
    for fa in fastas:
        n, bp = fasta_stats(fa) if os.path.getsize(fa) else (0, 0)
        print(f"{fa}\t{n}\t{bp}", file=log)
        if n:
            non_empty.append((fa, n, bp))
    big = [c for c in non_empty if c[2] >= min_bp]
    if big:
        print(f"{len(big)} of {len(fastas)} cells reach the minimum assembly "
              f"size ({min_bp:,} bp)", file=log)
        return big
    if non_empty:
        print(f"WARNING: no cell reaches the minimum assembly size "
              f"({min_bp:,} bp) -- choosing among all {len(non_empty)} "
              f"non-empty cells instead", file=log)
    return non_empty
