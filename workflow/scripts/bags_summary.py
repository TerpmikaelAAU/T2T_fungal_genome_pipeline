"""One annotation summary row per final assembly (see 11_annotation_bags.smk):
geneML gene count, BUSCO (protein mode), and antiSMASH regions by type --
the same numbers as BAGS' final report, plus a few more BGC classes.

BGC types come from each antiSMASH region GenBank's `/product=` qualifiers
on its `region` feature, so a hybrid region (e.g. NRPS + T1PKS) counts
towards every class it contains, and once in the total.
"""

import glob
import os
import re

snakemake = snakemake  # noqa: F821  (injected by Snakemake)

# antiSMASH product names -> the summary columns
CLASSES = {
    "NRPS":    lambda p: "NRPS" in p,
    "PKS":     lambda p: "PKS" in p,
    "terpene": lambda p: p == "terpene",
    "RiPP":    lambda p: "RiPP" in p or p in {"fungal-RiPP", "lanthipeptide", "RRE-containing"},
    "other":   None,  # anything not above
}


def region_products(gbk):
    """The /product= values of the `region` feature in one region GenBank."""
    products, in_region = [], False
    with open(gbk) as fh:
        for line in fh:
            if re.match(r"^ {5}\S", line):  # a new feature starts
                in_region = line.split()[0] == "region"
            elif in_region:
                m = re.search(r'/product="([^"]+)"', line)
                if m:
                    products.append(m.group(1))
    return products


genes = 0
with open(snakemake.input.gff) as fh:
    for line in fh:
        parts = line.split("\t")
        if len(parts) > 2 and parts[2] == "gene":
            genes += 1

busco = "NA"
with open(snakemake.input.busco) as fh:
    m = re.search(r"C:\S+", fh.read())
    if m:
        busco = m.group(0)

regions = sorted(glob.glob(os.path.join(snakemake.input.antismash, "*region*.gbk")))
counts = {c: 0 for c in CLASSES}
for gbk in regions:
    products = region_products(gbk)
    hit = False
    for cls, test in CLASSES.items():
        if test and any(test(p) for p in products):
            counts[cls] += 1
            hit = True
    if not hit:
        counts["other"] += 1

with open(snakemake.output.a, "w") as out:
    out.write("assembly\tgenes\tbusco_proteins\tbgc_regions\t"
              + "\t".join(f"{c}_regions" for c in CLASSES) + "\n")
    out.write(f"{snakemake.params.name}\t{genes}\t{busco}\t{len(regions)}\t"
              + "\t".join(str(counts[c]) for c in CLASSES) + "\n")
