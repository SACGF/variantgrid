#!/bin/bash
# Adds the GRCh38 reference allele sites gnomAD's T2T liftover drops to the T2T combined gnomAD file (#2788)
# @see gnomad_t2t_ref_swaps.py. Writes the 4.1.t2t_liftover_ref_swaps file next to the original, which it leaves alone.
#
# Usage: VEP_ANNOTATION_DIR=/data/annotation/VEP ANNOTATION_BASE_DIR=/data/annotation WORK_DIR=... \
#            gnomad_t2t_ref_swaps.sh [GRCh38 chromosome...]
# A test run names one chromosome, sets T2T_CONTIGS to its T2T contig and OUTPUT to somewhere else
# Needs bcftools, tabix, bgzip and a python3 with numpy + pysam (VariantGrid's virtualenv: set PYTHON)

set -e -o pipefail
export LC_ALL=C

for v in VEP_ANNOTATION_DIR ANNOTATION_BASE_DIR WORK_DIR; do
  if [[ -z ${!v} ]]; then
    echo "Please set '${v}'" >&2
    exit 1
  fi
done

GENERATE_DIR=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
VARIANTGRID_DIR=$(dirname $(dirname $(dirname ${GENERATE_DIR})))
PYTHON=${PYTHON:-python3}
THREADS=${THREADS:-4}

GNOMAD_GRCH38=${VEP_ANNOTATION_DIR}/annotation_data/GRCh38/gnomad4.1_GRCh38_contigs.vcf.gz
T2T_DIR=${VEP_ANNOTATION_DIR}/annotation_data/T2T-CHM13v2.0
GNOMAD_T2T=${T2T_DIR}/gnomad4.1.t2t_liftover_T2T-CHM13v2.0_combined_af.vcf.bgz
OUTPUT=${OUTPUT:-${T2T_DIR}/gnomad4.1.t2t_liftover_ref_swaps_T2T-CHM13v2.0_combined_af.vcf.bgz}
CHAIN=${ANNOTATION_BASE_DIR}/liftover/hg38ToHs1.over.chain.gz
GRCH38_FASTA=${ANNOTATION_BASE_DIR}/fasta/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
T2T_FASTA=${ANNOTATION_BASE_DIR}/fasta/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

if [[ $# -gt 0 ]]; then
  CHROMOSOMES="$@"
else
  CHROMOSOMES="$(seq 1 22) X Y"
fi

mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

# On a box shared with PostgreSQL / gunicorn / celery, each step first waits until MIN_FREE_GB (default 4) plus what
# it needs is available: ~1.5 GB for the sites step, ~0.5 GB for a swap step
MIN_FREE_GB=${MIN_FREE_GB:-4}
wait_for_memory() {
  local needed_gb=$1
  while (( $(awk '/MemAvailable/ {print int($2 / 1048576)}' /proc/meminfo) < MIN_FREE_GB + needed_gb )); do
    echo "Waiting for $((MIN_FREE_GB + needed_gb)) GB available memory..."
    sleep 60
  done
}

sites_for_chrom() {
  local chrom=$1
  if [[ ! -e sites_${chrom}.tsv ]]; then
    wait_for_memory 2
    ${PYTHON} -I ${GENERATE_DIR}/gnomad_t2t_ref_swaps.py sites --chain ${CHAIN} \
      --grch38-fasta ${GRCH38_FASTA} --t2t-fasta ${T2T_FASTA} \
      --grch38-chrom-mapping ${VARIANTGRID_DIR}/snpdb/genome/chrom_mapping_GRCh38.map \
      --t2t-chrom-mapping ${VARIANTGRID_DIR}/snpdb/genome/chrom_mapping_T2T-CHM13v2.0.map \
      --chrom ${chrom} --output sites_${chrom}.tmp
    mv sites_${chrom}.tmp sites_${chrom}.tsv
  fi
}

swaps_for_chrom() {
  local chrom=$1
  if [[ ! -e swaps_${chrom}.txt ]]; then
    wait_for_memory 1
    local contig=$(awk -v c=${chrom} '$1==c {print $2}' ${VARIANTGRID_DIR}/snpdb/genome/chrom_mapping_GRCh38.map)
    cut -f1,2 sites_${chrom}.tsv | uniq > targets_${chrom}.tsv
    bcftools view -r ${contig} -T targets_${chrom}.tsv ${GNOMAD_GRCH38} \
      | ${PYTHON} -I ${GENERATE_DIR}/gnomad_t2t_ref_swaps.py swap --sites sites_${chrom}.tsv \
          --existing-vcf ${GNOMAD_T2T} --output swaps_${chrom}.tmp 2> swaps_${chrom}.log
    mv swaps_${chrom}.tmp swaps_${chrom}.txt
  fi
}
export -f swaps_for_chrom wait_for_memory
export MIN_FREE_GB PYTHON GENERATE_DIR VARIANTGRID_DIR CHAIN GRCH38_FASTA T2T_FASTA GNOMAD_GRCH38 GNOMAD_T2T

# The sites step holds two chromosomes' sequences (~1GB+), so runs one at a time; the swap step is small
echo "Finding differing sites..."
for chrom in ${CHROMOSOMES}; do
  sites_for_chrom ${chrom}
done
echo "Finding swapped records..."
echo ${CHROMOSOMES} | tr ' ' '\n' | xargs -P ${THREADS} -I{} bash -c 'swaps_for_chrom {}'

echo "Sorting swapped records..."
SWAPS_VCF=swaps.vcf.gz
HEADER=header.vcf
(
  bcftools view -h ${GNOMAD_T2T} | grep -v "^#CHROM"
  # gnomad_data.py writes these without declaring them, which bcftools annotate (a backfill) needs
  echo '##INFO=<ID=AC_XY,Number=1,Type=Integer,Description="Alternate allele count for XY samples">'
  echo '##INFO=<ID=AN_XY,Number=1,Type=Integer,Description="Total number of alleles in XY samples">'
  echo '##INFO=<ID=AF_XY,Number=1,Type=Float,Description="Alternate allele frequency in XY samples">'
  echo "##gnomad_t2t_ref_swaps=Adds the GRCh38 reference allele at sites where CHM13 carries a GRCh38 ALT allele (and the site's other SNV alleles), lifted from $(basename ${GNOMAD_GRCH38}) by VariantGrid gnomad_t2t_ref_swaps.py (#2788)"
  bcftools view -h ${GNOMAD_T2T} | grep "^#CHROM"
) > ${HEADER}
(
  cat ${HEADER}
  for chrom in ${CHROMOSOMES}; do cat swaps_${chrom}.txt; done | sort -t$'\t' -k1,1 -k2,2n -k4,4 -k5,5 -u
) | bgzip -@ ${THREADS} > ${SWAPS_VCF}
tabix -f -p vcf ${SWAPS_VCF}

# A text merge: the T2T file has INFO fields its header doesn't define, which bcftools refuses to write back out
echo "Merging into the T2T liftover per contig..."
CONTIGS=${T2T_CONTIGS:-$(tabix -l ${GNOMAD_T2T})}  # T2T_CONTIGS: a test run's contigs
echo ${CONTIGS} | tr ' ' '\n' | xargs -P ${THREADS} -I{} \
  bash -c "[[ -e merged_{}.vcf.gz ]] || (sort -m -s -t\$'\\t' -k2,2n <(tabix ${GNOMAD_T2T} {}) <(tabix ${SWAPS_VCF} {}) | bgzip -@ 2 > merged_{}.tmp.vcf.gz && mv merged_{}.tmp.vcf.gz merged_{}.vcf.gz)"
# bgzip blocks concatenate - each chunk is removed once appended, so the merge needs one copy's disk, not two
bgzip < ${HEADER} > ${OUTPUT}.tmp
for contig in ${CONTIGS}; do
  cat merged_${contig}.vcf.gz >> ${OUTPUT}.tmp
  rm merged_${contig}.vcf.gz
done
mv ${OUTPUT}.tmp ${OUTPUT}
tabix -f -p vcf ${OUTPUT}
tabix -f --csi -p vcf ${OUTPUT}
echo "Wrote ${OUTPUT}"
