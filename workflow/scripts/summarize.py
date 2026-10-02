"""Build a single human-readable summary per sample.

Pulls together seqkit read stats, seqkit assembly stats, the BUSCO short
summary (before and after polishing, for polished samples), and a coverage
estimate derived from corrected bases / assembly size.
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
    if "hifiasm_telo" in path:
        return f"assembly {parent} +telo"
    if "hifiasm" in path:
        return f"assembly {parent}"
    if "contig" in path:
        return "selected assembly"
    # data/reads/ (FASTQ samples) or data/samtools/Fastq/ (BAM samples):
    # the input reads, minus the simplex parents of any duplex reads.
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

selector = snakemake.params.selector
selector_desc = {
    "lowest_contig": "fewest contigs",
    "highest_busco": "highest BUSCO completeness",
    "most_t2t": "most telomere-to-telomere contigs",
}.get(selector, selector)

out = []
out.append("=" * 70)
out.append(f"  SUMMARY: {snakemake.params.sample}  (selector: {selector} -- {selector_desc})")
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
out.append(f"ASSEMBLIES  (one per min_q/min_len grid cell; {selector_desc} is selected)")
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

def busco_section(title, busco_dir):
    """Append one BUSCO short summary to the report; return its C: percent
    (None if the summary can't be found or parsed)."""
    out.append("-" * 70)
    out.append(title)
    out.append("-" * 70)
    summary_file = find_busco(str(busco_dir))
    complete = None
    if summary_file:
        with open(summary_file) as fh:
            for line in fh:
                if re.search(r"C:|Complete|Fragmented|Missing|Total BUSCO|BUSCO skipped", line):
                    out.append("  " + line.strip())
                m = re.search(r"C:([\d.]+)%", line)
                if m and complete is None:
                    complete = float(m.group(1))
    else:
        out.append(f"  short_summary*.txt not found under {busco_dir}")
        out.append("  (check the BUSCO rule's -o path; see the run log)")
    out.append("")
    return complete


lineage = snakemake.params.lineage
if "busco_unpolished" in snakemake.input.keys():
    # pod5/bam samples: the same assembly before and after dorado polish.
    before = busco_section(f"BUSCO BEFORE POLISHING  ({lineage})",
                           snakemake.input.busco_unpolished)
    after = busco_section(f"BUSCO AFTER POLISHING  ({lineage})",
                          snakemake.input.busco)
    if before is not None and after is not None:
        out.append(f"Polishing changed BUSCO completeness: C {before:.1f}% -> "
                   f"{after:.1f}% ({after - before:+.1f})")
        out.append("")
else:
    busco_section(f"BUSCO  ({lineage})", snakemake.input.busco)

# --- telomeres (only with `telomere: enabled: true`) -----------------------
if "telomeres" in snakemake.input.keys():
    motif = {}
    with open(snakemake.input.motif) as fh:
        for line in fh:
            key, _, value = line.rstrip("\n").partition("\t")
            motif[key] = value
    contigs = []
    with open(snakemake.input.telomeres) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("contig\t"):
                continue
            contigs.append(line.rstrip("\n").split("\t"))
    ends = sum((c[4] == "yes") + (c[5] == "yes") for c in contigs)
    t2t = sum(c[6] == "yes" for c in contigs)
    out.append("-" * 70)
    out.append("TELOMERES")
    out.append("-" * 70)
    out.append(f"Motif             : {motif.get('motif', '?')}  "
               f"(5' end: {motif.get('motif_5prime', '?')})")
    out.append(f"hifiasm re-run    : {motif.get('rerun', '?')} -- {motif.get('reason', '')}")
    out.append(f"T2T contigs       : {t2t} of {len(contigs)}")
    out.append(f"Telomeric ends    : {ends} of {2 * len(contigs)}")
    out.append(f"  (per contig: {snakemake.input.telomeres})")
    out.append("")

# --- contamination (only with `contamination: enabled: true`) -------------
if "tiara" in snakemake.input.keys():
    lengths, name = {}, None
    with open(snakemake.input.final) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    classes = {}  # class -> [contigs, bp]
    prokaryotic = []
    with open(snakemake.input.tiara) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig, first, second = parts[:3]
            # organelle -> mitochondrion/plastid; everything else by stage 1
            cls = second if first == "organelle" and second != "n/a" else first
            bp = lengths.pop(contig, 0)
            classes.setdefault(cls, [0, 0])
            classes[cls][0] += 1
            classes[cls][1] += bp
            if first in ("bacteria", "archaea", "prokarya"):
                prokaryotic.append((bp, contig, first))
    if lengths:  # contigs Tiara skipped as too short
        classes["not classified (short)"] = [len(lengths), sum(lengths.values())]
    out.append("-" * 70)
    out.append("CONTAMINATION  (Tiara; contigs < "
               f"{snakemake.params.tiara_min_len:,} bp aren't classified)")
    out.append("-" * 70)
    out.append(f"{'class':<32}{'contigs':>12}{'bp':>16}")
    for cls, (n, bp) in sorted(classes.items(), key=lambda kv: -kv[1][1]):
        out.append(f"{cls:<32}{n:>12,}{bp:>16,}")
    if prokaryotic:
        prokaryotic.sort(reverse=True)
        out.append("")
        out.append(f"WARNING: {len(prokaryotic)} contig(s) classified as prokaryotic "
                   f"({sum(p[0] for p in prokaryotic):,} bp) -- possible contamination:")
        for bp, contig, cls in prokaryotic[:20]:
            out.append(f"  {contig:<30}{bp:>14,} bp  {cls}")
        if len(prokaryotic) > 20:
            out.append(f"  ... and {len(prokaryotic) - 20} more")
    else:
        out.append("")
        out.append("No contigs classified as bacteria or archaea.")
    out.append(f"  (per contig: {snakemake.input.tiara})")
    out.append("")

with open(snakemake.output.a, "w") as fh:
    fh.write("\n".join(out) + "\n")