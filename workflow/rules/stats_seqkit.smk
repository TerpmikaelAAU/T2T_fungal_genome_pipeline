configfile: "config/config.yaml"

rule seqkit_raw:
    input:
        "data/{input}.fastq"
    output:
        "data/seqkit/{input}/reads_raw.txt"
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["time"],
    conda:
       "seqkit.yml"
    shell:
        """
        seqkit stats {input} -Ta > {output}
        """

rule seqkit_porechopped:
    input:
        rules.porechop_abi.output.a
    output:
        "data/seqkit/{input}/reads_porechopped.txt"
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["time"],
    conda:
       "seqkit.yml"
    shell:
        """
        seqkit stats {input} -Ta > {output}
        """

rule seqkit_chopper_pre_dorado:
    input:
        rules.chopper_dorado.output.a
    output:
        "data/seqkit/{input}/reads_chopper_pre_dorado.txt"
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["time"],
    conda:
       "seqkit.yml"
    shell:
        """
        seqkit stats {input} -Ta > {output}
        """

rule seqkit_chopper_pre_flye:
    input:
        rules.chopper_flye.output.a
    output:
        "data/seqkit/{input}/reads_chopper_pre_flye.txt"
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["time"],
    conda:
       "seqkit.yml"
    shell:
        """
        seqkit stats {input} -Ta > {output}
        """

rule seqkit_chopper_ultralong:
    input:
        rules.chopper_ultralong.output.a
    output:
        "data/seqkit/{input}/reads_chopper_ultralong.txt"
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["time"],
    conda:
       "seqkit.yml"
    shell:
        """
        seqkit stats {input} -Ta > {output}
        """

rule seqkit_flye_assembly:
    input:
        rules.flye.output.b
    output:
        "data/seqkit/{input}/flye_assembly.txt"
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["time"],
    conda:
       "seqkit.yml"
    shell:
        """
        seqkit stats {input} -Ta > {output}
        """

rule seqkit_reads_after_correction:
    input:
        rules.dorado.output.a
    output:
        "data/seqkit/{input}/reads_after_correction.txt"
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["time"],
    conda:
       "seqkit.yml"
    shell:
        """
        seqkit stats {input} -Ta > {output}
        """

rule seqkit_reads_chopper_for_hifiasm:
    input:
        rules.chopper_hifiasm.output.a
    output:
        "data/seqkit/{input}/reads_chopper_for_hifiasm.txt"
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["time"],
    conda:
       "seqkit.yml"
    shell:
        """
        seqkit stats {input} -Ta > {output}
        """