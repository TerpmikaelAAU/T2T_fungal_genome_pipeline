configfile: "config/config.yaml"

rule contig_count:
    input:
        fasta_files = expand("data/hifiasm/{{input}}_{length}/{{input}}_{length}.fa", length=config["length"], input=config["input"])
    output:
        a = "data/contig/{input}_lowest_contig_file.fa"
    threads:
        5
    resources:
        mem_mb=resources["fga"]["mem_mb"],
        runtime=resources["fga"]["time"],
        partition="shared"
    shell:
      """
        mkdir -p $(dirname {output})
        lowest_contig_file=$(for fasta_file in {input.fasta_files}; do
            count=$(grep -c "^>" "$fasta_file" 2>/dev/null || echo -1)  # Count the number of contigs (lines starting with ">")
            echo "$count $fasta_file"
        done | sort -n | awk '$1 >= 0 {{print $2; exit}}')
        if [ -z "$lowest_contig_file" ]; then
            echo "No valid contig counts found." >&2
            exit 1
        fi
        cp "$lowest_contig_file" {output}
        """