# Last step for samples that can be polished (can_polish()): polish a
# winning assembly using the aligned reads (dorado_align). dorado picks the
# model from the basecaller model in the BAM's @RG header, and its more
# accurate move-table model only when the reads carry one (pod5 basecalled
# here, or a BAM/FASTQ basecalled with --emit-moves). GPU-only. Runs once per
# DISTINCT winner (see 6_05_polish_once.smk); polished_assembly then hands
# the result to every selector that picked it. Output is temp() --
# final_genome (9_stats.smk) copies it into
# results/{sample}/{selector}/{sample}_{selector}_final.fasta.
# GPU by default; config `dorado_polish.device: "cpu"` runs it on a CPU node
# instead (POLISH_DEVICE in the Snakefile). --threads matters there: without
# it dorado uses every core on the node, not the ones SLURM gave the job.
# --ignore-read-groups: a BAM merged from several runs (or a run that was
# restarted) has one read group per run, which polish otherwise rejects;
# all of them must still share one basecalling model.
rule dorado_polish:
    input:
        dorado = dorado_bin,
        a = rules.dorado_align.output.a,
        bai = rules.dorado_align.output.bai,  # dorado polish needs the index
        b = unique_winner,
    output:
        a = temp("data/dorado_polish_unique/{input}_{selector}.fasta")
    threads:
        POLISH_THREADS
    resources:
        mem_mb=resources["dorado_polish"]["mem_mb"],
        runtime=resources["dorado_polish"]["runtime"],
        **({"gres": POLISH_GRES} if POLISH_GRES else {}),
    params:
        infer_threads = lambda w, threads: threads if POLISH_DEVICE == "cpu" else 2,
    log:
        "logs/dorado_polish/{input}_{selector}.log"
    shell:
        """
        "{input.dorado}" polish --batchsize 8 --device {POLISH_DEVICE} \
            --threads {threads} --infer-threads {params.infer_threads} \
            --ignore-read-groups {input.a} {input.b} > {output.a} 2> {log}
        """