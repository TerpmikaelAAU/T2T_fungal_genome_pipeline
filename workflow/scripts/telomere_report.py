"""Per-contig telomere table for one final assembly:
results/{sample}/{selector}/telomeres.tsv (see 9_1_telomere_report.smk).
summarize.py turns it into the TELOMERES section of summary.txt."""

import sys

snakemake = snakemake  # noqa: F821  (injected by Snakemake)
sys.path.insert(0, snakemake.scriptdir)
from telomeres import contig_ends, read_search  # noqa: E402

motif = "?"
with open(snakemake.input.motif) as fh:
    for line in fh:
        key, _, value = line.rstrip("\n").partition("\t")
        if key == "motif":
            motif = value

rows = contig_ends(read_search(snakemake.input.search), snakemake.params.min_repeats)
rows.sort(key=lambda r: -r["length"])
with open(snakemake.output.a, "w") as out:
    out.write(f"# motif: {motif}; an end is telomeric with >= "
              f"{snakemake.params.min_repeats} repeats in its terminal "
              f"{snakemake.params.window} bp window\n")
    out.write("contig\tlength\tstart_repeats\tend_repeats\ttelomere_start\ttelomere_end\tt2t\n")
    for r in rows:
        out.write(f"{r['contig']}\t{r['length']}\t{r['start_repeats']}\t{r['end_repeats']}\t"
                  f"{'yes' if r['telomere_start'] else 'no'}\t"
                  f"{'yes' if r['telomere_end'] else 'no'}\t"
                  f"{'yes' if r['telomere_start'] and r['telomere_end'] else 'no'}\n")
