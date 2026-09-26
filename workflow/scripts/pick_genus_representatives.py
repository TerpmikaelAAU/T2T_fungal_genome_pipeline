"""Pick ONE NCBI assembly per genus to serve as the phylogeny's reference set
(see 0_4_ncbi_genus_representatives.smk).

Genus comes from NCBI taxonomy (the lineage lookup done by
ncbi_genome_summary). Only when an organism's taxid is missing from that
lookup does it fall back to the first word of the organism name -- and
names that clearly aren't a genus ("uncultured fungus", "Ascomycota sp.")
are skipped rather than guessed at.

Within a genus, the representative is the first after sorting by:
  1. NCBI "reference genome" designation first
  2. assembly level: Complete Genome > Chromosome > Scaffold > Contig
  3. contig N50, highest first
  4. release date, newest first

Output: one row per genus, tip label first (that's also the genome's file
name after download, and so its label in the tree), accession second.
"""

import json
import re

snakemake = snakemake  # noqa: F821  (injected by Snakemake)

LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
RANKS = ["phylum", "class", "order", "family", "genus"]

# Placeholder / above-genus names that turn up as the first word of an
# organism name ("Fungi sp.", "Ascomycota sp. XYZ", "Mucorales sp.", ...).
NOT_A_GENUS = {"Fungi", "Eukaryota", "Dikarya", "Opisthokonta", "Incertae"}
HIGHER_RANK_SUFFIXES = ("mycota", "mycotina", "mycetes", "mycetidae",
                        "ales", "aceae", "ineae", "oideae")

COLUMNS = ["label", "accession", "genus", "organism", "phylum", "class",
           "order", "family", "assembly_level", "refseq_category",
           "contig_n50", "release_date"]


def load_lineages(path):
    """taxid -> {rank: name} from `datasets summary taxonomy --as-json-lines`."""
    lineages = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            record = json.loads(line)
            tax = record.get("taxonomy", record)
            taxid = tax.get("tax_id", tax.get("taxId"))
            if taxid is None:
                continue
            classification = tax.get("classification") or {}
            lineages[int(taxid)] = {
                rank: (classification.get(rank) or {}).get("name", "")
                for rank in RANKS
            }
    return lineages


def genus_from_name(name):
    """Fallback: first word of the organism name, if it looks like a genus."""
    words = name.replace("[", "").replace("]", "").split()
    if not words:
        return ""
    word = words[0]
    if (not word.isalpha() or not word[0].isupper() or word in NOT_A_GENUS
            or word.endswith(HIGHER_RANK_SUFFIXES)):
        return ""
    return word


def tip_label(organism, accession):
    """Newick/filename-safe: Genus_species_strain_GCA_000000000.1"""
    name = re.sub(r"[^A-Za-z0-9]+", "_", organism).strip("_")[:60].rstrip("_")
    return f"{name}_{accession}" if name else accession


def rank_key(row):
    date = int(re.sub(r"\D", "", row["release_date"]) or 0)
    return (
        0 if row["refseq_category"] == "reference genome" else 1,
        LEVEL_RANK.get(row["assembly_level"], len(LEVEL_RANK)),
        -int(row["contig_n50"] or 0),
        -date,
    )


def main():
    lineages = load_lineages(snakemake.input.taxonomy)

    best = {}  # genus -> row
    n_assemblies = n_no_genus = n_not_current = 0
    with open(snakemake.input.genomes) as fh:
        for line in fh:
            if not line.strip():
                continue
            n_assemblies += 1
            report = json.loads(line)
            organism = report.get("organism", {})
            info = report.get("assemblyInfo", {})
            stats = report.get("assemblyStats", {})

            status = info.get("assemblyStatus", "current")
            if status != "current":
                n_not_current += 1
                continue

            name = organism.get("organismName", "")
            taxid = organism.get("taxId")
            lineage = lineages.get(int(taxid)) if taxid is not None else None
            if lineage is not None:
                genus = lineage["genus"]  # "" if NCBI has no genus for it
            else:
                lineage = {rank: "" for rank in RANKS}
                genus = genus_from_name(name)
            if not genus:
                n_no_genus += 1
                continue

            accession = report["accession"]
            row = {
                "label": tip_label(name, accession),
                "accession": accession,
                "genus": genus,
                "organism": name,
                "phylum": lineage["phylum"],
                "class": lineage["class"],
                "order": lineage["order"],
                "family": lineage["family"],
                "assembly_level": info.get("assemblyLevel", ""),
                "refseq_category": info.get("refseqCategory", ""),
                "contig_n50": str(stats.get("contigN50", "")),
                "release_date": info.get("releaseDate", ""),
            }
            if genus not in best or rank_key(row) < rank_key(best[genus]):
                best[genus] = row

    if not best:
        raise ValueError(
            f"No usable assemblies in {snakemake.input.genomes} -- check "
            "`phylogeny.taxon` in config.yaml and "
            "logs/phylogeny/ncbi_genome_summary.log."
        )

    rows = sorted(best.values(),
                  key=lambda r: [r[rank] for rank in RANKS[:-1]] + [r["genus"]])
    with open(snakemake.output.table, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            # Tabs/newlines inside NCBI names would break the TSV.
            out.write("\t".join(re.sub(r"\s", " ", row[c]) for c in COLUMNS) + "\n")

    with open(snakemake.log[0], "w") as log:
        print(f"{n_assemblies} assemblies read; "
              f"{n_not_current} not current, {n_no_genus} without a genus "
              f"(skipped); {len(rows)} genera -> {snakemake.output.table}",
              file=log)
        if not lineages:
            print("WARNING: no taxonomy lineages available -- genus taken "
                  "from organism names and phylum..family left blank",
                  file=log)


main()
