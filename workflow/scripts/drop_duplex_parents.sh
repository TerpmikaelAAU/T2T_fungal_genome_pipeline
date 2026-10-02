#!/usr/bin/env bash
# Plain-text copy of a FASTQ sample's reads (decompressed once if .gz), with
# the simplex parents of duplex reads removed -- used by prepare_fastq
# (0_3_decompress.smk).
#
# dorado duplex writes every duplex read (named "<template id>;<complement
# id>") AND its two simplex parent reads, so the same molecule would enter
# hifiasm three times. The parents are dropped by name, keeping the duplex
# read and every simplex read without a duplex partner. Simplex-only data has
# no "A;B" names, so it passes through unchanged.
#
# Usage: drop_duplex_parents.sh <in.fastq[.gz]> <out.fastq> <threads>
set -euo pipefail
in=$1
out=$2
threads=${3:-1}
tmp="$out.tmp"
parents="$out.parents"
trap 'rm -f "$tmp" "$parents"' EXIT

if [[ $in == *.gz ]]; then
    if command -v pigz >/dev/null 2>&1; then
        pigz -dc -p "$threads" "$in" > "$tmp"
    else
        zcat "$in" > "$tmp"
    fi
    src=$tmp
else
    src=$in
fi

# Pass 1: count reads and duplex reads; write both parent ids of each duplex
# read to $parents.
read -r total duplex < <(awk -v parents="$parents" '
    NR % 4 == 1 {
        id = substr($1, 2)
        if (index(id, ";")) {
            duplex++
            n = split(id, p, ";")
            for (i = 1; i <= n; i++) print p[i] > parents
        }
    }
    END { close(parents); print NR / 4, duplex + 0 }' "$src")

if [ "$duplex" -eq 0 ]; then
    echo "$total reads, no duplex reads (no read names like A;B) -- all reads kept"
    if [ "$src" = "$tmp" ]; then mv "$tmp" "$out"; else cp "$src" "$out"; fi
else
    # Pass 2: copy every record whose name isn't a duplex parent.
    awk 'NR == FNR { drop[$1] = 1; next }
         FNR % 4 == 1 { keep = !(substr($1, 2) in drop) }
         keep' "$parents" "$src" > "$out"
    kept=$(( $(wc -l < "$out") / 4 ))
    echo "$total reads, $duplex duplex: dropped $(( total - kept )) simplex parent reads of duplex reads, $kept reads kept"
fi
