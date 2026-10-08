# Last step for samples that can be polished (can_polish()): polish a
# winning assembly using the aligned reads (dorado_align). dorado picks the
# model from the basecaller model in the BAM's @RG header, and its more
# accurate move-table model only when the reads carry one (pod5 basecalled
# here, or a BAM/FASTQ basecalled with --emit-moves). GPU-only. Runs once per
# DISTINCT winner (see 6_05_polish_once.smk); polished_assembly then hands
# the result to every selector that picked it. Output is temp() --
# final_genome (9_stats.smk) copies it into
# results/{sample}/{selector}/{sample}_{selector}_final.fasta.
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
        16
    resources:
        mem_mb=resources["dorado_polish"]["mem_mb"],
        runtime=resources["dorado_polish"]["runtime"],
        gres=GPU_GRES,
    shell:
        """
        "{input.dorado}" polish --batchsize 8 --device cuda:all --ignore-read-groups {input.a} {input.b} > {output.a}
        
        """