# ============================================================================
#  Polish each distinct winner only once
# ============================================================================
# The selectors (lowest_contig, highest_busco, most_t2t -- see Snakefile
# SELECTORS) often pick the SAME grid cell. Aligning and polishing it once
# per selector would repeat the most expensive GPU step for an identical
# result, so for samples that can be polished (can_polish()):
#   1. polish_groups (checkpoint) -- compares the selectors' winners byte for
#      byte and keeps one copy of each distinct assembly, named after the
#      first selector (in SELECTORS order) that picked it. The mapping is
#      written to results/{sample}/polish_groups.tsv.
#   2. dorado_align / dorado_polish -- run once per distinct assembly
#      (6_1_dorado_align.smk, 6_3_dorado_polish.smk).
#   3. polished_assembly -- gives every selector its polished assembly under
#      data/dorado_polish/{sample}_{selector}.fasta, as before, so BUSCO,
#      final_genome and everything downstream are unchanged.

checkpoint polish_groups:
    input:
        winners = lambda w: [f"data/contig/{w.input}_{sel}_file.fa" for sel in SELECTORS],
    output:
        # Not temp(): dorado_align reads from here, and it's small (one
        # unpolished assembly per distinct winner).
        dir    = directory("data/polish_groups/{input}"),
        groups = "results/{input}/polish_groups.tsv",
    params:
        selectors = SELECTORS,
    threads:
        1
    resources:
        mem_mb  = 2000,
        runtime = 30,
    run:
        import hashlib, shutil
        os.makedirs(output.dir, exist_ok=True)
        polished_as = {}  # md5 -> selector whose copy gets polished
        with open(output.groups, "w") as fh:
            fh.write("selector\tpolished_as\tmd5\n")
            for sel, fa in zip(params.selectors, input.winners):
                with open(fa, "rb") as f:
                    md5 = hashlib.md5(f.read()).hexdigest()
                if md5 not in polished_as:
                    polished_as[md5] = sel
                    shutil.copy(fa, os.path.join(output.dir, f"{sel}.fa"))
                fh.write(f"{sel}\t{polished_as[md5]}\t{md5}\n")


def polished_as(wildcards):
    """The selector whose copy of this selector's winner gets polished --
    itself, unless an earlier selector picked the identical assembly."""
    groups = checkpoints.polish_groups.get(input=wildcards.input).output.groups
    with open(groups) as fh:
        next(fh)
        for line in fh:
            sel, canon, _ = line.rstrip("\n").split("\t")
            if sel == wildcards.selector:
                return canon
    raise WorkflowError(f"{wildcards.selector} missing from {groups}")


def unique_winner(wildcards):
    """dorado_align's draft: this distinct assembly's copy in polish_groups."""
    d = checkpoints.polish_groups.get(input=wildcards.input).output.dir
    return os.path.join(d, f"{wildcards.selector}.fa")


rule polished_assembly:
    input:
        a = lambda w: f"data/dorado_polish_unique/{w.input}_{polished_as(w)}.fasta",
    output:
        a = temp("data/dorado_polish/{input}_{selector}.fasta"),
    threads:
        1
    resources:
        mem_mb  = 2000,
        runtime = 30,
    shell:
        """
        cp {input.a} {output.a}
        """
