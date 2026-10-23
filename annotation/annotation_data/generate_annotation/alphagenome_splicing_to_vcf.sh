#!/bin/bash
# AlphaGenome merged splicing scores (#2433) as a VCF for VEP --custom, which only matches alleles for VCF.
#
# Input: combined_alphagenome_splicing_snvs.tsv.gz (+ .tbi), extracted from combined_splicing_snvs_tabix.zip
# ("AlphaGenome SNV merged splicing scores", https://deepmind.google.com/science/alphagenome/downloads - Google
# sign-in, non-commercial use only). The zip is stored, not compressed, so extracting it is a plain copy. Rows are
# hg38 SNVs on chr1-22, X and Y: '#CHROM POS REF ALT alphagenome_splicing'.
#
# Output: one record per row with the score as INFO/score, contigs without 'chr' to match our other GRCh38 custom
# VCFs, bgzipped + tabix indexed. Chromosomes are converted in parallel, then the bgzip parts concatenated in order.
#
# Usage: alphagenome_splicing_to_vcf.sh <combined_alphagenome_splicing_snvs.tsv.gz> <output.vcf.gz> [jobs]

set -e -o pipefail

INPUT=$1
OUTPUT=$2
JOBS=${3:-3}

if [[ -z "${INPUT}" || -z "${OUTPUT}" ]]; then
  echo "Usage: $0 <combined_alphagenome_splicing_snvs.tsv.gz> <output.vcf.gz> [jobs]" >&2
  exit 1
fi
if [[ ! -e "${INPUT}.tbi" ]]; then
  echo "Missing tabix index ${INPUT}.tbi" >&2
  exit 1
fi

PARTS_DIR=$(mktemp -d "$(dirname "${OUTPUT}")/alphagenome_splicing_parts.XXXXXX")
trap 'rm -rf "${PARTS_DIR}"' EXIT

CONTIGS=$(tabix -l "${INPUT}")

{
  echo "##fileformat=VCFv4.2"
  echo "##source=alphagenome_splicing_to_vcf.sh $(basename "${INPUT}")"
  for contig in ${CONTIGS}; do
    echo "##contig=<ID=${contig#chr}>"
  done
  echo '##INFO=<ID=score,Number=1,Type=Float,Description="AlphaGenome merged splicing score (non-commercial use only)">'
  printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
} | bgzip -c > "${PARTS_DIR}/000.vcf.gz"

convert_contig() {
  local order=$1 contig=$2
  tabix "${INPUT}" "${contig}" \
    | mawk -F '\t' -v OFS='\t' '{ sub(/^chr/, "", $1); print $1, $2, ".", $3, $4, ".", ".", "score=" $5 }' \
    | bgzip -c > "${PARTS_DIR}/${order}.vcf.gz"
  echo "${contig} done" >&2
}
export -f convert_contig
export INPUT PARTS_DIR

order=0
for contig in ${CONTIGS}; do
  order=$((order + 1))
  printf "%03d %s\n" ${order} "${contig}"
done | xargs -P "${JOBS}" -L 1 bash -c 'convert_contig "$0" "$1"'

# bgzip files concatenate into a valid bgzip file; the zero-padded names sort in contig order
cat "${PARTS_DIR}"/*.vcf.gz > "${OUTPUT}"
tabix -p vcf "${OUTPUT}"
echo "Wrote ${OUTPUT}: $(tabix -l "${OUTPUT}" | wc -l) contigs"
