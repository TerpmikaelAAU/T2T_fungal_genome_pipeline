"""Build one sample's results/<sample>/report.html (rule report, 9_2_report.smk).

Collects the numbers behind the page -- what the read filters kept, each
grid cell's assembly, BUSCO and hifiasm k-mer histogram, and which cell each
selector picked -- into one JSON object and drops it into
report_template.html, which draws the charts in the browser.

Forgiving like summarize.py: a missing hifiasm log (cells assembled before
hifiasm got a log file) or an unreadable BUSCO summary leaves a gap on the
page instead of failing the run.
"""

import bisect
import datetime
import glob
import gzip
import json
import math
import os
import re

snakemake = snakemake  # noqa: F821  (injected by Snakemake)
P = snakemake.params


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


# --- reads: a length x mean-quality histogram, in bases ---------------------
# Every cutoff the grid uses is a bin edge, so the page can add up bins and
# get exactly what each cell keeps (chopper keeps length >= min_len and mean
# Q >= min_q, and seqkit computes mean Q the same way chopper does).
def length_edges(max_len, cutoffs):
    top = max(max_len + 1, 1001)
    edges = {round(10 ** (2 + i / 20)) for i in range(int(20 * (math.log10(top) - 2)) + 2)}
    return sorted(edges | {int(c) for c in cutoffs})


def q_edges(max_q, cutoffs):
    return sorted({i / 2 for i in range(int(2 * max_q) + 2)} | {float(c) for c in cutoffs})


def read_table(path):
    """data/report/<sample>/{reads,corrected}.tsv.gz: one 'length<TAB>mean Q'
    line per read (just the length for the corrected reads). Streamed, not
    held in memory: a run can have tens of millions of reads."""
    if not path:
        return
    n = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            # A malformed line (e.g. seqkit's columns changing) stops the job
            # here rather than drawing empty charts.
            f = line.rstrip("\n").split("\t")
            yield int(f[0]), (float(f[1]) if len(f) > 1 else 0.0)
            n += 1
    if not n:
        raise SystemExit(f"{path} lists no reads")


def bin_index(edges, value):
    return min(max(bisect.bisect_right(edges, value) - 1, 0), len(edges) - 2)


reads_path = snakemake.input.reads_hist
corrected_path = snakemake.input.corrected_hist if P.dorado_correct else None
len_cutoffs = list(P.grid["min_len"]) + ([P.dorado_len] if P.dorado_correct else [])
q_cutoffs = list(P.grid["min_q"]) + ([P.dorado_q] if P.dorado_correct else [])
# Two passes: the first finds how far the bins must reach.
max_len, max_q = 1000, 1.0
for length, q in read_table(reads_path):
    max_len, max_q = max(max_len, length), max(max_q, q)
for length, _ in read_table(corrected_path):
    max_len = max(max_len, length)
le = length_edges(max_len, len_cutoffs)
qe = q_edges(max_q, q_cutoffs)
hist = [[0] * (len(qe) - 1) for _ in range(len(le) - 1)]
for length, q in read_table(reads_path):
    hist[bin_index(le, length)][bin_index(qe, q)] += length
corrected_hist = [0] * (len(le) - 1)
for length, _ in read_table(corrected_path):
    corrected_hist[bin_index(le, length)] += length


# --- read stats: the raw reads, and what hifiasm got per cell ----------------
stats = {row.get("file", ""): row for row in read_tsv(snakemake.input.read_stats)}
raw = next((row for f, row in stats.items()
            if f.startswith(("data/reads/", "data/samtools/Fastq/"))), {})


# --- hifiasm k-mer histogram -------------------------------------------------
# hifiasm prints one histogram after k-mer counting and one after each
# round of read correction; the last one is what it assembles and purges
# duplicates from, and its homozygous coverage matches that histogram's peak.
HIST = re.compile(r"\[M::ha_hist_line\]\s+(\d+|rest):\s*\**>?\s*(\d+)\s*$")
SIDE = re.compile(r"\[M::ha_analyze_count\] (left|right): (?:count\[(\d+)\]|none)")
HOM = re.compile(r"homozygous read coverage threshold: (\d+)")


def kmer_histogram(log):
    try:
        with open(log, errors="replace") as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None
    blocks, current, side, hom = [], None, {}, None
    for line in lines:
        m = HIST.search(line)
        if m:
            if current is None:
                current, side = {}, {}
                blocks.append((current, side))
            current[m[1]] = int(m[2])
            continue
        current = None
        s = SIDE.search(line)
        if s and blocks:
            blocks[-1][1][s[1]] = int(s[2]) if s[2] else None
        h = HOM.search(line)
        if h:
            hom = int(h[1])
    if not blocks:
        return None
    counts, side = blocks[-1]
    depths = sorted(int(k) for k in counts if k != "rest")
    if not depths:
        return None
    series = [counts.get(str(d), 0) for d in range(1, depths[-1] + 1)]
    peak = hom
    if not peak:
        # Highest bar past the error spike: after the first local minimum.
        start = next((i for i in range(1, len(series) - 1)
                      if series[i] <= series[i + 1]), 0)
        peak = max(range(start, len(series)), key=series.__getitem__) + 1
    # left = a peak below the main one (heterozygous); right = one above it.
    second = side.get("left") or side.get("right")
    return {"counts": series, "rest": counts.get("rest", 0), "peak": peak, "het": second}


cells = []
for c in P.cells:
    row = stats.get(c["reads"], {})
    cells.append({
        "q": c["q"], "l": c["l"],
        "reads": as_int(row.get("num_seqs")), "bases": as_int(row.get("sum_len")),
        "n50": as_int(row.get("N50")), "capped": bool(P.capped),
        "kmer": kmer_histogram(c["hifiasm_log"]),
    })


# --- assemblies: size, contigs, BUSCO, T2T per grid cell ---------------------
asm_stats = {row.get("file", ""): row for row in read_tsv(snakemake.input.asm_stats)}
BUSCO = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)")


def parse_busco(text):
    """The 'C:..%[S:..%,D:..%],F:..%,M:..%,n:..' line; None when BUSCO didn't
    run (the n:0 placeholder for an assembly too small to search)."""
    m = BUSCO.search(text or "")
    if not m or int(m[6]) == 0:
        return None
    return {"S": float(m[2]), "D": float(m[3]), "F": float(m[4]), "M": float(m[5]), "n": int(m[6])}


grid_busco = {}
with open(snakemake.input.busco_grid) as fh:
    for line in fh:
        d, _, text = line.rstrip("\n").partition("\t")
        grid_busco[d] = parse_busco(text)


def t2t_counts(log):
    """pick_most_t2t.py's per-cell table: fasta, t2t_contigs, ends, contigs."""
    out = {}
    try:
        with open(log) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == 4 and f[1].isdigit():
                    out[f[0]] = int(f[1])
    except OSError:
        pass
    return out


t2t = t2t_counts(P.select_logs["most_t2t"]) if "most_t2t" in P.select_logs else {}
assemblies = []
for a in P.assemblies:
    row = asm_stats.get(a["fasta"], {})
    assemblies.append({
        "pass": a["pass"], "q": a["q"], "l": a["l"],
        "asm_len": as_int(row.get("sum_len")), "contigs": as_int(row.get("num_seqs")),
        "busco": grid_busco.get(a["busco_dir"]),
        "t2t": t2t.get(a["fasta"]) if t2t else None,
    })


# --- which cell each selector picked -----------------------------------------
by_fasta = {a["fasta"]: a for a in P.assemblies}
winners = {}
for sel, log in P.select_logs.items():
    try:
        with open(log) as fh:
            picked = [ln.split("picked:", 1)[1].strip() for ln in fh if ln.startswith("picked:")]
    except OSError:
        picked = []
    if picked and picked[-1] in by_fasta:
        a = by_fasta[picked[-1]]
        winners[sel] = {"pass": a["pass"], "q": a["q"], "l": a["l"]}


def busco_of(busco_dir):
    hits = glob.glob(os.path.join(busco_dir, "**", "short_summary*.txt"), recursive=True)
    if not hits:
        return None
    with open(hits[0]) as fh:
        return parse_busco(fh.read())


finals = []
for sel in P.selectors:
    if P.polished:
        finals.append({"selector": sel, "stage": "before polishing",
                       "busco": busco_of(snakemake.input.busco_unpolished[P.selectors.index(sel)])})
        finals.append({"selector": sel, "stage": "polished",
                       "busco": busco_of(snakemake.input.busco_final[P.selectors.index(sel)])})
    else:
        finals.append({"selector": sel, "stage": "",
                       "busco": busco_of(snakemake.input.busco_final[P.selectors.index(sel)])})


# --- genome size, for coverage -----------------------------------------------
# The highest-BUSCO pick: the assembly most likely to hold the whole genome.
genome_size, size_from = 0, ""
for sel in ("highest_busco", "lowest_contig"):
    w = winners.get(sel)
    if w:
        a = next(a for a in assemblies if (a["pass"], a["q"], a["l"]) == (w["pass"], w["q"], w["l"]))
        if a["asm_len"]:
            genome_size = a["asm_len"]
            size_from = "size of the " + sel.replace("_", "-").replace("highest-busco", "highest-BUSCO") + " assembly"
            break
if not genome_size:
    genome_size = max([a["asm_len"] for a in assemblies] + [1])
    size_from = "largest grid assembly"

data = {
    "sample": P.sample,
    "generated": datetime.date.today().isoformat(),
    "genome_size": genome_size, "genome_size_from": size_from,
    "dorado_correct": bool(P.dorado_correct), "dorado_q": P.dorado_q, "dorado_len": P.dorado_len,
    "capped": bool(P.capped), "polished": bool(P.polished), "busco_lineage": P.lineage,
    "grid": {"min_q": list(P.grid["min_q"]), "min_len": list(P.grid["min_len"])},
    "reads_total": {"bases": as_int(raw.get("sum_len")), "reads": as_int(raw.get("num_seqs")),
                    "n50": as_int(raw.get("N50"))},
    "reads": {"len_edges": le, "q_edges": qe, "hist": hist, "corrected": corrected_hist},
    "cells": cells, "assemblies": assemblies, "winners": winners, "finals": finals,
}

with open(os.path.join(snakemake.scriptdir, "report_template.html")) as fh:
    page = fh.read()
# "</" can't appear inside the <script> block; JSON allows escaping the slash.
blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
with open(snakemake.output.a, "w") as fh:
    fh.write(page.replace("/*REPORT_DATA*/null", blob, 1))
