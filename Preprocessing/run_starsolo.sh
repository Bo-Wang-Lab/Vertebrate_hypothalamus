#!/bin/bash
#
# run_starsolo.sh — loop STARsolo over Danio_rerio_hypo* folders
#
# Auto-detects the *_1 (CB+UMI) and *_2 (cDNA) fastqs in each folder and
# comma-joins them in matching order. Works whether a folder has one pair
# or two (lanes of the same library).
#
# Usage:  ./run_starsolo.sh            # all folders
#         ./run_starsolo.sh Danio_rerio_hypo9   # one folder

set -euo pipefail

GENOME=DR_genome_ncbi/
WHITELIST=737K-august-2016.txt
THREADS=8

# folders: either the ones passed as args, or all matching the glob
if [[ $# -gt 0 ]]; then
    folders=("$@")
else
    folders=(Danio_rerio_hypo*/)
fi

for d in "${folders[@]}"; do
    d="${d%/}"                                   # strip trailing slash
    [[ -d "$d" ]] || { echo "[skip] $d: not a directory"; continue; }

    # --- collect reads, sorted so _1 and _2 stay in matching order --------
    mapfile -t r1s < <(ls "$d"/*_1.fastq.gz 2>/dev/null | sort)
    mapfile -t r2s < <(ls "$d"/*_2.fastq.gz 2>/dev/null | sort)

    if [[ ${#r1s[@]} -eq 0 ]]; then
        echo "[skip] $d: no *_1.fastq.gz found"; continue
    fi
    if [[ ${#r1s[@]} -ne ${#r2s[@]} ]]; then
        echo "[SKIP] $d: ${#r1s[@]} R1 vs ${#r2s[@]} R2 — mismatched, not running"
        continue
    fi

    # comma-join (no spaces) for STAR
    R1=$(IFS=,; echo "${r1s[*]}")
    R2=$(IFS=,; echo "${r2s[*]}")

    # --- skip if already done ---------------------------------------------
    if [[ -s "$d/Solo.out/Gene/Summary.csv" ]]; then
        echo "[done] $d: Solo.out already present, skipping"; continue
    fi

    echo "======================================================"
    echo "[run ] $d  (${#r1s[@]} pair(s))"
    printf '       R2: %s\n       R1: %s\n' "$R2" "$R1"
    echo "[time] start $(date +%T)"

    STAR --genomeDir "$GENOME" \
        --readFilesIn "$R2" "$R1" \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
        --soloUMIlen 10 \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \
        --soloFeatures Gene \
        --soloMultiMappers EM \
        --outSAMtype None \
        --runThreadN "$THREADS" \
        --outFilterMultimapNmax 20 \
        --readFilesCommand zcat \
        --outFileNamePrefix "$d"/

    echo "[time] done  $(date +%T)"
done

# --- collect key metrics across all folders -------------------------------
echo
echo "=== Gene summary across folders ==="
for d in "${folders[@]}"; do
    d="${d%/}"
    s="$d/Solo.out/Gene/Summary.csv"
    [[ -s "$s" ]] || continue
    cells=$(awk -F, '/Estimated Number of Cells/{print $2}' "$s")
    med=$(awk -F, '/Median UMI per Cell/{print $2}' "$s")
    gene=$(awk -F, '/Reads Mapped to Gene: Unique Gene/{print $2}' "$s")
    printf "%-24s cells=%-8s medUMI=%-8s uniqGene=%s\n" "$d" "$cells" "$med" "$gene"
done
