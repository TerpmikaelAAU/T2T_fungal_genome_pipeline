# "fga" = FASTA-from-GFA: hifiasm's primary-contig graph (prefix.p_ctg.gfa)
# encodes each contig as an 'S' (segment) line; this pulls just the
# name/sequence pairs into a plain FASTA for everything downstream.
rule fga:
    input:
        a = rules.hifiasm.output.a,
    output:
        a =  temp("data/hifiasm/{input}_q{minq}_l{minlen}/{input}_q{minq}_l{minlen}.fa")
    threads:
        12
    resources:
        mem_mb=resources["fga"]["mem_mb"],
        runtime=resources["fga"]["runtime"],
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.a} > {output.a}

        """