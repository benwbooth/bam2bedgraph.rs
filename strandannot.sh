#!/usr/bin/env bash
# Given an GFF/GTF file and a tab-delimited genome sizes file,
# creates a strand annotation bam file for use with bam2bedgraph
# --autostrand option. Requires bedtools and samtools to be installed.

set -eu -o pipefail
if [[ $# < 2 ]]; then
  echo "Usage: $0 gff-file genome-file" >&2
  exit 1
fi

bedtools sort -i "$1" >"${1%.*}.sorted.$(basename "${1##*.}")"
bedtools merge -c 7 -o distinct -s -i "${1%.*}.sorted.$(basename "${1##*.}")" |perl -F\\t -lane 'print join "\t",$F[0],"strandannot","merged_$.",$F[1]+1,$F[2],".",$F[3],".","gene_id \"merged_$.\"; transcript_id \"merged_$.\";"' >"${1%.*}.sorted.merged.$(basename "${1##*.}")"
bedtools subtract -S -a "${1%.*}.sorted.merged.$(basename "${1##*.}")" -b "${1%.*}.sorted.merged.$(basename "${1##*.}")" >"${1%.*}.autostrand.gtf"
#bedtools subtract -S -a "${1%.*}.sorted.merged.$(basename "${1##*.}")" -b "${1%.*}.sorted.merged.$(basename "${1##*.}")" >"${1%.*}.sorted.merged.subtracted.$(basename "${1##*.}")"
#bedtools bedtobam -i "${1%.*}.sorted.merged.subtracted.$(basename "${1##*.}")" -g "$2" >"${1%.*}.autostrand.bam"
#samtools index "${1%.*}.autostrand.bam"
