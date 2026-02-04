#!/bin/bash

## Before running, process fastqs with nf-core methylseq
## /home/rzhang/nextflow run /path/ro/methylseq/ --input /path/to/sample_sheets.csv  --outdir /path/to/output/ --genome mm10 -profile docker,gpu --aligner bwameth   -resume


set -euo pipefail

# Path to mm10 chromosome sizes
MM10_CHROM_SIZES="/data/rzhang/mm10_em_control/genome.chrom.sizes"

for bedgraph in ./BedMethyl/*.bedGraph; do
    [[ -f "$bedgraph" ]] || continue

    base="${bedgraph%.bedGraph}"

    echo "Processing $bedgraph ..."

    # Remove header if any
    awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $4, $5, $6}' "$bedgraph" > "${base}.tmp"

    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4/100}' "${base}.tmp" > "${base}.methylfrac.bedGraph"

    # Create methylated count bedGraph
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' "${base}.tmp" > "${base}.methyl_count.bedGraph"

    # Create coverage bedGraph (methyl + unmethyl)
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5+$6}' "${base}.tmp" > "${base}.coverage.bedGraph"

    # Sort bedGraphs
    sort -k1,1 -k2,2n "${base}.methylfrac.bedGraph" > "${base}.methylfrac.sorted.bedGraph"
    sort -k1,1 -k2,2n "${base}.methyl_count.bedGraph" > "${base}.methyl_count.sorted.bedGraph"
    sort -k1,1 -k2,2n "${base}.coverage.bedGraph" > "${base}.coverage.sorted.bedGraph"

    # Convert to bigWig
    bedGraphToBigWig "${base}.methyl_count.sorted.bedGraph" "$MM10_CHROM_SIZES" "${base}.methyl_count.bw"
    bedGraphToBigWig "${base}.coverage.sorted.bedGraph" "$MM10_CHROM_SIZES" "${base}.coverage.bw"
    bedGraphToBigWig "${base}.methylfrac.sorted.bedGraph" "$MM10_CHROM_SIZES" "${base}.methylfrac.bw"
    

    # Clean up intermediate files
    rm -f "${base}.tmp" "${base}.methyl_count.bedGraph" "${base}.coverage.bedGraph" \
          "${base}.methyl_count.sorted.bedGraph" "${base}.coverage.sorted.bedGraph" \
          "${base}.methylfrac.sorted.bedGraph" "${base}.methylfrac.bedGraph"

    echo "Generated: ${base}.methyl_count.bw and ${base}.coverage.bw"
done

echo "✅ All done."
