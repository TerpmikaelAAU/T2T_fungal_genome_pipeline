#!/usr/bin/env bash
# Download the reads of one or more SRA/ENA/DDBJ runs from ENA as a single
# .fastq.gz -- used by download_sra (0_5_download_sra.smk).
#
# For each run, ENA's file report gives the sequencing platform, the FASTQ
# URL(s) and their md5s. Each file is downloaded over https, checked against
# its md5, and the files are joined in the order given (concatenated gzip
# files are one valid gzip file).
#
# Usage: download_sra.sh <out.fastq.gz> <run accession> [run accession ...]
set -euo pipefail
out=$1
shift
part="$out.part"
trap 'rm -f "$part" "$part".*' EXIT

api="https://www.ebi.ac.uk/ena/portal/api/filereport"
fields="run_accession,instrument_platform,instrument_model,read_count,base_count,fastq_ftp,fastq_md5"
files=()
for acc in "$@"; do
    report=$(curl -fsS --retry 5 --retry-delay 30 \
        "$api?accession=$acc&result=read_run&fields=$fields&format=tsv")
    # Pick columns by header name rather than trusting their order.
    row=$(awk -F'\t' -v fields="$fields" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        NR == 2 { n = split(fields, f, ",")
                  for (i = 1; i <= n; i++) printf "%s%s", $col[f[i]], (i < n ? "|" : "\n") }' \
        <<< "$report")
    if [ -z "$row" ]; then
        echo "ERROR: ENA has no run $acc. Check the accession; a run made public only in the last few days may not have reached ENA yet."
        exit 1
    fi
    # "|", not a tab: read merges runs of tabs, so an empty field would shift.
    IFS='|' read -r run platform model reads bases ftp md5 <<< "$row"
    echo "$run: $platform $model, $reads reads, $bases bases"
    if [ "$platform" != OXFORD_NANOPORE ]; then
        echo "ERROR: $run was sequenced on $platform, not Oxford Nanopore; this pipeline assembles Nanopore reads only."
        exit 1
    fi
    if [ -z "$ftp" ]; then
        echo "ERROR: ENA has no FASTQ for $run yet. Download it with fasterq-dump (sra-tools) and add it as a fastq sample instead."
        exit 1
    fi
    IFS=';' read -ra urls <<< "$ftp"
    IFS=';' read -ra sums <<< "$md5"
    for i in "${!urls[@]}"; do
        f="$part.${#files[@]}"
        echo "Downloading https://${urls[$i]}"
        curl -fsSL --retry 5 --retry-delay 30 "https://${urls[$i]}" -o "$f"
        got=$(md5sum "$f" | cut -d' ' -f1)
        if [ "$got" != "${sums[$i]}" ]; then
            echo "ERROR: md5 of ${urls[$i]} is $got, ENA lists ${sums[$i]}: the download is incomplete or corrupt. Rerun to try again."
            exit 1
        fi
        files+=("$f")
    done
done

if [ "${#files[@]}" -eq 1 ]; then
    mv "${files[0]}" "$part"
else
    cat "${files[@]}" > "$part"
fi
mv "$part" "$out"
echo "$(du -h "$out" | cut -f1) written to $out from ${#files[@]} file(s)"
