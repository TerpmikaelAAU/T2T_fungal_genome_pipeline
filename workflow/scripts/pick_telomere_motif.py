"""Pick a sample's telomere motif from the `tidk explore` candidates and
decide whether hifiasm should be re-run with it (see
5_1_telomere_motif.smk).

Each of the top candidates from `tidk explore` was searched with
`tidk search` on the same assembly. The best motif is the one that makes
the most contig ends telomeric (ties: higher explore count, then shorter).
It "covers" a fraction of the telomeres: its telomeric ends out of every
end any candidate found telomeric. hifiasm is re-run with it
(--telo-m) only if that fraction and its number of ends both clear the
config thresholds -- i.e. one motif clearly explains the telomeres.

If explore found nothing usable, the configured fallback motif is used
for counting telomeres downstream, and hifiasm is not re-run.
"""

import os
import sys

snakemake = snakemake  # noqa: F821  (injected by Snakemake)
sys.path.insert(0, snakemake.scriptdir)
from telomeres import contig_ends, read_search, revcomp, tally  # noqa: E402

p = snakemake.params
cand_dir = snakemake.input.candidates

candidates = []
with open(os.path.join(cand_dir, "candidates.tsv")) as fh:
    next(fh)
    for line in fh:
        rank, motif, count = line.rstrip("\n").split("\t")[:3]
        rows = contig_ends(
            read_search(os.path.join(cand_dir, f"rank{rank}_telomeric_repeat_windows.tsv")),
            p.min_repeats)
        ends, t2t, _ = tally(rows)
        telomeric = {(r["contig"], side) for r in rows
                     for side in ("start", "end") if r[f"telomere_{side}"]}
        candidates.append({"rank": int(rank), "motif": motif, "explore_count": int(count),
                           "ends": ends, "t2t": t2t, "telomeric": telomeric, "rows": rows})

all_ends = set().union(*(c["telomeric"] for c in candidates)) if candidates else set()
best = max(candidates, key=lambda c: (c["ends"], c["explore_count"], -len(c["motif"])),
           default=None)

if best and best["ends"] > 0:
    motif = best["motif"]
    fraction = best["ends"] / len(all_ends)
    # hifiasm wants the motif as it reads at a chromosome's 5' (left) end:
    # whichever orientation dominates the first window of telomeric starts.
    fwd = sum(r["start_forward"] for r in best["rows"] if r["telomere_start"])
    rev = sum(r["start_reverse"] for r in best["rows"] if r["telomere_start"])
    motif_5prime = motif if fwd >= rev else revcomp(motif)
    rerun = (p.rerun_enabled and best["ends"] >= p.min_ends
             and fraction >= p.min_fraction)
    reason = (f"best motif covers {fraction:.0%} of {len(all_ends)} telomeric ends "
              f"found by any candidate ({best['ends']} ends; need >= "
              f"{p.min_fraction:.0%} and >= {p.min_ends} ends)")
else:
    motif = p.fallback
    # Written the conventional way (3' end, e.g. TTAGGG), so the 5' end
    # reads as its reverse complement. Only informational: no re-run here.
    motif_5prime = revcomp(p.fallback)
    fraction = 0.0
    rerun = False
    reason = f"no candidate motif found telomeres; using fallback {p.fallback}"
if not p.rerun_enabled:
    reason += "; re-run disabled in config"

with open(snakemake.output.motif, "w") as out:
    for key, value in [("motif", motif), ("motif_5prime", motif_5prime),
                       ("fraction", f"{fraction:.3f}"),
                       ("rerun", "yes" if rerun else "no"), ("reason", reason)]:
        out.write(f"{key}\t{value}\n")

with open(snakemake.output.table, "w") as out:
    out.write("rank\tmotif\texplore_count\ttelomeric_ends\tt2t_contigs\tchosen\n")
    for c in candidates:
        out.write(f"{c['rank']}\t{c['motif']}\t{c['explore_count']}\t{c['ends']}\t"
                  f"{c['t2t']}\t{'yes' if c is best and best['ends'] > 0 else ''}\n")
    out.write(f"# motif used: {motif} (5' end: {motif_5prime}); "
              f"hifiasm re-run: {'yes' if rerun else 'no'} -- {reason}\n")
