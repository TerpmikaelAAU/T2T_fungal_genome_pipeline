rule hifiasm:
    input:
        unpack(hifiasm_inputs)
    output:
        a = temp("data/hifiasm/{input}_q{minq}_l{minlen}/prefix.p_ctg.gfa"),
        b = temp(["data/hifiasm/{input}_q{minq}_l{minlen}/prefix.a_ctg.gfa",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.a_ctg.lowQ.bed",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.a_ctg.noseq.gfa",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.p_ctg.lowQ.bed",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.p_ctg.noseq.gfa",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.p_utg.gfa",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.p_utg.lowQ.bed",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.p_utg.noseq.gfa",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.r_utg.gfa",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.r_utg.lowQ.bed",
        "data/hifiasm/{input}_q{minq}_l{minlen}/prefix.r_utg.noseq.gfa",
        ])

    threads:
        12
    resources:
        mem_mb=scaled_mem(6.0, 64000),
        runtime=scaled_time(0.3, 720),
    params:
        ul_flag    = lambda w, input: f"--ul {input.b}" if config["ultralong"]["enabled"] else "",
        min_bytes  = int(config.get("hifiasm_min_input_mb", 1)) * 1_000_000,
        prefix     = lambda w: f"data/hifiasm/{w.input}_q{w.minq}_l{w.minlen}/prefix",
    conda:
        "../envs/hifiasm.yml"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        size=$(stat -c%s {input.a} 2>/dev/null || echo 0)
        if [ "$size" -lt {params.min_bytes} ]; then
            # A near-empty grid point SIGILLs hifiasm instead of failing cleanly
            # (see Handoff.md §5.2). Write empty placeholders so contig_count's
            # own filter (count > 0) skips this cell instead of retrying it.
            echo "WARNING: {input.a} is ${{size}} bytes (< {params.min_bytes}); skipping hifiasm, writing empty placeholders" >&2
            : > {output.a}
            for f in {output.b}; do : > "$f"; done
        else
            hifiasm --ont -t {threads} --primary {params.ul_flag} -o "{params.prefix}" {input.a}
        fi
        """
