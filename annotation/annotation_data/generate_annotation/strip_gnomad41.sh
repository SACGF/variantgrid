#!/bin/bash
set -euo pipefail

#
# Generate one per-chromosome script that strips a gnomAD v4.1 joint VCF
# down to the INFO fields we need. Each script calls strip_gnomad41_one.sh
# for its chromosome and also renames contigs (NC_000001.11 -> 1 etc) using
# the supplied bcftools chrom map.
#
# Behaviour:
#   - No PARTITION argument: write the per-chrom scripts and stop. The
#     user runs them however they like (gnu parallel, plain bash loop,
#     manual sbatch, etc).
#   - PARTITION supplied: write the scripts with a #SBATCH --partition=...
#     directive AND submit them via sbatch.
#
# Usage:
#   bash strip_gnomad41.sh INPUT_DIR OUTPUT_DIR CHROM_MAP [PARTITION]
#
# CHROM_MAP is a bcftools --rename-chrs TSV, e.g.
#   snpdb/genome/chrom_mapping_GRCh38.map
#

USAGE="Usage: $0 INPUT_DIR OUTPUT_DIR CHROM_MAP [PARTITION]"
INPUT_DIR="${1:?$USAGE}"
OUTPUT_DIR="${2:?$USAGE}"
CHROM_MAP="${3:?$USAGE}"
PARTITION="${4:-}"

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

# PARTITION presence drives both the #SBATCH directive and whether we
# auto-submit. Without it we still emit valid SLURM scripts (just minus
# the partition line, so a later manual `sbatch` would fall back to the
# cluster default — or the user can edit the file).
if [ -n "$PARTITION" ]; then
    PARTITION_DIRECTIVE="#SBATCH --partition=${PARTITION}"
    SUBMIT=1
else
    PARTITION_DIRECTIVE=""
    SUBMIT=0
fi

SCRIPTS=()
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
${PARTITION_DIRECTIVE}
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
export PATH="${EXTRA_PATH}:\${PATH}"
echo "Processing chr${CHR} - \$(date)"

RENAME_FILE="${RENAME_FILE}" NON_PAR_HEADER_FILE="${NON_PAR_HEADER_FILE}" bash "${WORKER}" "${INPUT_VCF}" "${OUTPUT_VCF}" "${CHROM_MAP}"

echo "Done chr${CHR} - \$(date)"
SLURM
    chmod +x "$SLURM_SCRIPT"
    SCRIPTS+=("$SLURM_SCRIPT")

    if [ "$SUBMIT" = "1" ]; then
        JID=$(sbatch --parsable "$SLURM_SCRIPT")
        echo "  Submitted chr${CHR} -> job $JID"
        JOBIDS+=($JID)
    else
        echo "  Wrote $SLURM_SCRIPT"
    fi
done

echo ""
if [ "$SUBMIT" = "1" ]; then
    echo "Submitted ${#JOBIDS[@]} jobs to partition '${PARTITION}'"
    echo "Monitor with: squeue -u $USER"
else
    echo "Wrote ${#SCRIPTS[@]} scripts to $OUTPUT_DIR"
    echo ""
    echo "Run them however you like, e.g.:"
    echo "  sbatch:        re-run this script with a PARTITION argument"
    echo "  gnu parallel:  parallel --jobs 4 bash {} ::: $OUTPUT_DIR/strip_chr*.slurm"
    echo "  serial loop:   for s in $OUTPUT_DIR/strip_chr*.slurm; do bash \"\$s\"; done"
fi
