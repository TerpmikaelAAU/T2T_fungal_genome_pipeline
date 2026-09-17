# Picks the winning assembly across a sample's (min_q, min_len) grid: the
# fewest-contigs candidate, on the assumption that's the most contiguous
# assembly. CAVEAT for repeat-rich genomes (e.g. oomycetes): fewest contigs
# can mean COLLAPSED repeats rather than a genuinely better assembly --
# cross-check the winner's total length in assembly_stats.tsv against the
# expected genome size before trusting this pick.
def contig_count_candidates(wildcards):
    grid = filter_grid(wildcards.input)
    return expand("data/hifiasm/{input}_q{minq}_l{minlen}/{input}_q{minq}_l{minlen}.fa",
                  input=wildcards.input,
                  minq=grid["min_q"],
                  minlen=grid["min_len"])

rule contig_count:
    input:
        fasta_files = contig_count_candidates
    output:
        a = temp("data/contig/{input}_lowest_contig_file.fa")
    threads:
        5
    resources:
        mem_mb=resources["fga"]["mem_mb"],
        runtime=resources["fga"]["runtime"],
    shell:
      """
        mkdir -p $(dirname {output})
        n_candidates=$(echo {input.fasta_files} | wc -w)
        lowest_contig_file=$(for fasta_file in {input.fasta_files}; do
            count=$(grep -c "^>" "$fasta_file" 2>/dev/null || echo -1)  # Count the number of contigs (lines starting with ">")
            echo "$count $fasta_file"
        done | sort -n | awk '$1 > 0 {{print $2; exit}}')
        if [ -z "$lowest_contig_file" ]; then
            echo "No valid contig counts found -- every grid point produced an empty or missing assembly." >&2
            exit 1
        fi
        if [ "$n_candidates" -eq 1 ]; then
            echo "Note: only 1 grid cell configured for this sample -- picked it directly, nothing to compare against." >&2
        fi
        cp "$lowest_contig_file" {output}
        """
