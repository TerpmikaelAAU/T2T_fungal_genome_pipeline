rule fga:
    input:
        a = rules.hifiasm.output.a,
    output: 
        a =  temp("data/hifiasm/{input}_{length}/{input}_{length}.fa")
    threads:
        12
    resources:
        mem_mb=resources["fga"]["mem_mb"],
        runtime=resources["fga"]["runtime"],
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.a} > {output.a}

        """