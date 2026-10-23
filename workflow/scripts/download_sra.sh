#!/usr/bin/env bash
# Download the reads of one or more SRA/ENA/DDBJ runs from ENA as a single
# .fastq.gz -- used by download_sra (0_5_download_sra.smk).
#
# For each run, ENA's file report gives the sequencing platform, the FASTQ
# URL(s), their sizes and md5s. Each file is downloaded, checked against its
# md5, and the files are joined in the order given (concatenated gzip files
# are one valid gzip file).
#
# ENA's file server is slow (~1-2 MB/s) and drops long transfers, so each
# file is downloaded to <out>.partial.<file name> and resumed where it
# stopped: after a dropped connection, a stall, or an error such as a 403,
# the download carries on from the bytes already on disk. It gives up only
# after MAX_IDLE attempts in a row that add nothing; the job's walltime
# bounds the total. The partial files are kept when the script fails, so a
# retry or a later run also resumes, and removed once the output is
# complete.
#
# Usage: download_sra.sh <out.fastq.gz> <run accession> [run accession ...]
set -euo pipefail
out=$1
shift
joined="$out.part"
trap 'rm -f "$joined"' EXIT

MAX_IDLE=10        # attempts in a row without progress before giving up
HTTPS_FIRST=3      # failed https attempts in a row before trying ftp
WAIT=${DOWNLOAD_SRA_WAIT:-60}   # seconds between attempts

size_of() { if [ -e "$1" ]; then stat -c %s "$1"; else echo 0; fi; }

# fetch <url without scheme> <file> <expected bytes>
fetch() {
    local url=$1 f=$2 want=$3 idle=0 https_fails=0 tries=0 restarted=0
    local scheme=https have before rc
    while true; do
        have=$(size_of "$f")
        if [ "$have" -gt "$want" ]; then
            # Once only: if a download from zero overshoots too, ENA's
            # listed size is wrong, and starting again would loop forever.
            if [ "$restarted" -eq 1 ]; then
                rm -f "$f"
                echo "ERROR: ENA serves more than the $want bytes it lists for $url. Its listing is out of date; try again in a few days, or download the run with fasterq-dump and add it as a fastq sample."
                return 1
            fi
            echo "  $f is larger than the $want bytes ENA lists; starting it again"
            rm -f "$f"
            have=0
            restarted=1
        fi
        [ "$have" -eq "$want" ] && return 0
        if [ "$idle" -ge "$MAX_IDLE" ]; then
            echo "ERROR: no progress in $MAX_IDLE attempts in a row ($have of $want bytes). Rerun later to resume."
            return 1
        fi
        tries=$(( tries + 1 ))
        echo "  attempt $tries ($scheme): resuming at $have of $want bytes"
        before=$have
        # --speed-limit/--speed-time: drop a connection that has stalled
        # (< 10 kB/s for 2 min) instead of waiting out the walltime on it.
        rc=0
        curl -fsSL -C - --connect-timeout 60 --speed-limit 10000 --speed-time 120 \
            "$scheme://$url" -o "$f" || rc=$?
        [ "$rc" -ne 0 ] && echo "  curl exit $rc"
        if [ "$(size_of "$f")" -gt "$before" ]; then
            idle=0
            [ "$scheme" = https ] && https_fails=0
            continue
        fi
        # curl succeeding without adding anything means the server has no
        # bytes past ours (a 416): its file is smaller than ENA lists.
        if [ "$rc" -eq 0 ]; then
            echo "ERROR: the server has only $(size_of "$f") of the $want bytes ENA lists for $url. Its listing is out of date; try again in a few days, or download the run with fasterq-dump and add it as a fastq sample."
            return 1
        fi
        idle=$(( idle + 1 ))
        # Stay on whichever protocol works; once https has failed
        # HTTPS_FIRST times in a row, switch protocol after every failure.
        [ "$scheme" = https ] && https_fails=$(( https_fails + 1 ))
        if [ "$https_fails" -ge "$HTTPS_FIRST" ]; then
            if [ "$scheme" = https ]; then scheme=ftp; else scheme=https; fi
        fi
        sleep "$WAIT"
    done
}

api="https://www.ebi.ac.uk/ena/portal/api/filereport"
fields="run_accession,instrument_platform,instrument_model,read_count,base_count,fastq_ftp,fastq_bytes,fastq_md5"
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
    IFS='|' read -r run platform model reads bases ftp bytes md5 <<< "$row"
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
    IFS=';' read -ra sizes <<< "$bytes"
    IFS=';' read -ra sums <<< "$md5"
    if [ "${#sizes[@]}" -ne "${#urls[@]}" ] || [ "${#sums[@]}" -ne "${#urls[@]}" ] \
            || ! [[ "${sizes[*]}" =~ ^[0-9]+( [0-9]+)*$ ]]; then
        echo "ERROR: ENA lists ${#urls[@]} file(s) for $run but sizes '$bytes' and md5s '$md5'; can't check the download. Try again later."
        exit 1
    fi
    for i in "${!urls[@]}"; do
        f="$out.partial.${urls[$i]##*/}"
        echo "Downloading ${urls[$i]} (${sizes[$i]} bytes)"
        fetch "${urls[$i]}" "$f" "${sizes[$i]}"
        got=$(md5sum "$f" | cut -d' ' -f1)
        if [ "$got" != "${sums[$i]}" ]; then
            echo "  md5 is $got, ENA lists ${sums[$i]}; downloading it again from the start"
            rm -f "$f"
            fetch "${urls[$i]}" "$f" "${sizes[$i]}"
            got=$(md5sum "$f" | cut -d' ' -f1)
            if [ "$got" != "${sums[$i]}" ]; then
                rm -f "$f"
                echo "ERROR: md5 of ${urls[$i]} is $got twice, ENA lists ${sums[$i]}."
                exit 1
            fi
        fi
        echo "  md5 OK"
        files+=("$f")
    done
done

if [ "${#files[@]}" -eq 1 ]; then
    mv "${files[0]}" "$out"
else
    cat "${files[@]}" > "$joined"
    mv "$joined" "$out"
fi
# Also clears partials left over from an earlier, different accession list.
rm -f "$out".partial.*
echo "$(du -h "$out" | cut -f1) written to $out from ${#files[@]} file(s)"
