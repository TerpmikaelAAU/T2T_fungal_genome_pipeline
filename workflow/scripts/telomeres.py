"""Shared helpers for reading `tidk search` output (imported by
pick_telomere_motif.py, pick_most_t2t.py, telomere_report.py).

`tidk search` counts a motif (forward) and its reverse complement (reverse)
in fixed windows along each contig. A contig END counts as telomeric when
its terminal window holds at least `min_repeats` copies, in either
orientation; a contig with both ends telomeric is telomere-to-telomere
(T2T).
"""

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_search(path):
    """tidk search TSV -> {contig: [(window_end, forward, reverse), ...]},
    windows in file (= positional) order. Empty/missing -> {}."""
    contigs = {}
    try:
        fh = open(path)
    except FileNotFoundError:
        return contigs
    with fh:
        next(fh, None)  # header: id window forward_repeat_number ...
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            contigs.setdefault(parts[0], []).append(
                (int(parts[1]), int(parts[2]), int(parts[3])))
    return contigs


def contig_ends(contigs, min_repeats):
    """One row per contig: length, repeats in the first/last window, and
    whether each end is telomeric."""
    rows = []
    for contig, windows in contigs.items():
        first, last = windows[0], windows[-1]
        start = first[1] + first[2]
        end = last[1] + last[2]
        rows.append({
            "contig": contig,
            "length": last[0],  # the last window ends at the contig's end
            "start_repeats": start,
            "end_repeats": end,
            "telomere_start": start >= min_repeats,
            "telomere_end": end >= min_repeats,
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
