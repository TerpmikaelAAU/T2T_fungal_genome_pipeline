"""Shared helpers for reading `tidk search` output (imported by
pick_telomere_motif.py, pick_most_t2t.py, telomere_report.py).

`tidk search` counts a motif (forward) and its reverse complement (reverse)
in fixed windows along each contig. A contig END counts as telomeric when
its terminal window, both orientations together, holds at least `min_bp`
bp of motif AND at least `min_enrichment` times the contig's median
interior window. A contig with both ends telomeric is
telomere-to-telomere (T2T).

The bp floor stops short motifs from winning on copy count alone. The
enrichment test stops motifs that are common everywhere: in an AT-rich
genome, AAAAT reaches 60 bp per 2 kb window inside contigs as often as at
their ends, while a real telomere repeat is absent inside.
"""

import math
import statistics

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_search(path):
    """tidk search TSV -> ({contig: [(window_end, forward, reverse), ...]},
    motif), windows in file (= positional) order. Empty/missing -> ({}, "")."""
    contigs, motif = {}, ""
    try:
        fh = open(path)
    except FileNotFoundError:
        return contigs, motif
    with fh:
        next(fh, None)  # header: id window forward_repeat_number ... telomeric_repeat
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            contigs.setdefault(parts[0], []).append(
                (int(parts[1]), int(parts[2]), int(parts[3])))
            if len(parts) > 4:
                motif = parts[4]
    return contigs, motif


def contig_ends(search, min_bp, min_enrichment):
    """One row per contig of a read_search() result: length, repeats in the
    first/last window, and whether each end is telomeric (module docstring)."""
    contigs, motif = search
    min_copies = math.ceil(min_bp / max(len(motif), 1))
    rows = []
    for contig, windows in contigs.items():
        first, last = windows[0], windows[-1]
        start = first[1] + first[2]
        end = last[1] + last[2]
        interior = [f + r for _, f, r in windows[1:-1]]
        floor = max(min_copies,
                    min_enrichment * max(statistics.median(interior) if interior else 0, 1))
        rows.append({
            "contig": contig,
            "length": last[0],  # the last window ends at the contig's end
            "start_repeats": start,
            "end_repeats": end,
            "telomere_start": start >= floor,
            "telomere_end": end >= floor,
            # 5'-end orientation evidence, for hifiasm --telo-m
            "start_forward": first[1],
            "start_reverse": first[2],
        })
    return rows


def tally(rows):
    """(telomeric ends, T2T contigs, contigs)"""
    ends = sum(r["telomere_start"] + r["telomere_end"] for r in rows)
    t2t = sum(r["telomere_start"] and r["telomere_end"] for r in rows)
    return ends, t2t, len(rows)
