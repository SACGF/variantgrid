#!/bin/bash
set -euo pipefail

#
# Submit one SLURM job per chromosome to strip gnomAD v4.1 joint VCFs.
# Each job calls strip_gnomad41_one.sh for its chromosome, which also renames
# contigs (NC_000001.11 -> 1 etc) using the supplied bcftools chrom map.
#
# Usage:
#   bash strip_gnomad41.sh INPUT_DIR OUTPUT_DIR PARTITION CHROM_MAP [--dry-run]
#
# CHROM_MAP is a bcftools --rename-chrs TSV, e.g.
#   snpdb/genome/chrom_mapping_GRCh38.map
#

INPUT_DIR="${1:?Usage: $0 INPUT_DIR OUTPUT_DIR PARTITION CHROM_MAP [--dry-run]}"
OUTPUT_DIR="${2:?Usage: $0 INPUT_DIR OUTPUT_DIR PARTITION CHROM_MAP [--dry-run]}"
PARTITION="${3:?Usage: $0 INPUT_DIR OUTPUT_DIR PARTITION CHROM_MAP [--dry-run]}"
CHROM_MAP="${4:?Usage: $0 INPUT_DIR OUTPUT_DIR PARTITION CHROM_MAP [--dry-run]}"
DRY_RUN="${5:-}"

if [ ! -f "$CHROM_MAP" ]; then
    echo "CHROM_MAP not found: $CHROM_MAP" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER="${SCRIPT_DIR}/strip_gnomad41_one.sh"
# shellcheck source=gnomad41_fields.sh
source "${SCRIPT_DIR}/gnomad41_fields.sh"

# HPC tool paths (directories containing bcftools, tabix, bgzip etc.)
TOOL_PATHS=(
    /hpcfs/groups/phoenix-hpc-sacgf/tools/bcftools/current
    /hpcfs/groups/phoenix-hpc-sacgf/tools/tabix-0.2.6
)

CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

mkdir -p "$OUTPUT_DIR"

# Write rename file + non_par header file once; all per-chrom jobs share them.
RENAME_FILE="$OUTPUT_DIR/rename_annots.txt"
write_rename_file "$RENAME_FILE"
echo "Wrote $RENAME_FILE"

NON_PAR_HEADER_FILE="$OUTPUT_DIR/non_par.hdr"
write_non_par_header_file "$NON_PAR_HEADER_FILE"
echo "Wrote $NON_PAR_HEADER_FILE"

EXTRA_PATH=$(IFS=:; echo "${TOOL_PATHS[*]}")

JOBIDS=()
for CHR in $CHROMOSOMES; do
    INPUT_VCF="$INPUT_DIR/gnomad.joint.v4.1.sites.chr${CHR}.vcf.bgz"
    OUTPUT_VCF="$OUTPUT_DIR/gnomad4.1_chr${CHR}.stripped.vcf.gz"
    SLURM_SCRIPT="$OUTPUT_DIR/strip_chr${CHR}.slurm"

    cat > "$SLURM_SCRIPT" <<SLURM
#!/bin/bash
#SBATCH --job-name=gnomad41_chr${CHR}
#SBATCH --output=${OUTPUT_DIR}/gnomad41_chr${CHR}_%j.out
#SBATCH --error=${OUTPUT_DIR}/gnomad41_chr${CHR}_%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=${PARTITION}
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
export PATH="${EXTRA_PATH}:\${PATH}"
echo "Processing chr${CHR} - \$(date)"

RENAME_FILE="${RENAME_FILE}" NON_PAR_HEADER_FILE="${NON_PAR_HEADER_FILE}" bash "${WORKER}" "${INPUT_VCF}" "${OUTPUT_VCF}" "${CHROM_MAP}"

echo "Done chr${CHR} - \$(date)"
SLURM

    if [ "$DRY_RUN" = "--dry-run" ]; then
        echo "  [dry-run] Would submit $SLURM_SCRIPT"
    else
        JID=$(sbatch --parsable "$SLURM_SCRIPT")
        echo "  Submitted chr${CHR} -> job $JID"
        JOBIDS+=($JID)
    fi
done

echo ""
if [ "$DRY_RUN" = "--dry-run" ]; then
    echo "Dry run complete. 24 scripts written to $OUTPUT_DIR"
    echo "Review a script, then re-run without --dry-run to submit."
else
    echo "Submitted ${#JOBIDS[@]} jobs"
    echo "Monitor with: squeue -u $USER"
fi
