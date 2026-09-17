# Last step for pod5/bam entry points: polish the winning assembly (for
# whichever selector is being built -- see Snakefile SELECTORS) using the
# aligned reads (dorado_align) and their move table. GPU-only. Output is
# temp() -- final_genome (9_stats.smk) copies it straight into
# results/{sample}/{selector}/{sample}_{selector}_final.fasta, so this copy
# is redundant once that exists.
rule dorado_polish:
    input:
        dorado = dorado_bin,
        a = rules.dorado_align.output.a,
        b = "data/contig/{input}_{selector}_file.fa",
    output:
        a = temp("data/dorado_polish/{input}_{selector}.fasta")
    threads:
        16
    resources:
        mem_mb=resources["dorado_polish"]["mem_mb"],
        runtime=resources["dorado_polish"]["runtime"],
        gres=GPU_GRES,
    shell:
        """
        "{input.dorado}" polish --batchsize 8 --device cuda:all {input.a} {input.b} > {output.a}
        
        """