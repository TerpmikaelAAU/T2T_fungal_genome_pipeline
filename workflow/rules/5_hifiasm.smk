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
        32
    resources:
        # factor 4.0, ~20-25% headroom above observed peak RSS/input on the
        # p_infestans_88069 2x2 grid (actual ratio 2.4-3.3x, worst case at
        # the highest-coverage cells); was 6.0 and ran at 40-55% memory
        # utilization (efficiency_report_c170c0c2-*.csv).
        mem_mb=scaled_mem(4.0, 64000),
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
            # hifiasm's overlap/index caches (prefix.*.bin) aren't declared
            # Snakemake outputs, so temp() can't clean them -- and they're
            # large (100s of GB per grid cell). Nothing downstream reads
            # them once the .gfa files above exist, so delete them here.
            rm -f {params.prefix}.*.bin
        fi
        """
